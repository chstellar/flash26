import torch
import argparse
import os
import sys
import yaml
from tqdm import tqdm
import json
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp
from collections import defaultdict
import gc

sys.path.append(os.environ.get("SAFARI_PATH", "."))
from src.models.sequence.long_conv_lm import ConvLMHeadModel
from src.dataloaders.datasets.hg38_char_tokenizer import CharacterTokenizer

try:
    from tokenizers import Tokenizer
except:
    pass


class FastHG38Encoder:
    def __init__(
        self, model_cfg, ckpt_path, max_seq_len, nlayer, use_multiple_gpus=True
    ):
        self.max_seq_len = max_seq_len
        self.nlayer = nlayer
        self.use_multiple_gpus = use_multiple_gpus
        self.model, self.tokenizer = self.load_model(
            model_cfg, ckpt_path, max_seq_len, nlayer
        )

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.num_gpus = torch.cuda.device_count()

        if self.use_multiple_gpus and self.num_gpus > 1:
            print(f"Using {self.num_gpus} GPUs for inference")
            self.model = torch.nn.DataParallel(self.model)
        else:
            print(f"Using single GPU/CPU for inference")

        self.model = self.model.to(self.device)
        self.model.eval()  # Set to eval mode

    def tokenize_sequences_parallel(self, sequences, num_workers=4):
        """Parallel tokenization of sequences"""

        def tokenize_batch(seq_batch):
            tokenized = []
            for seq in seq_batch:
                if isinstance(self.tokenizer, Tokenizer):
                    tokens = self.tokenizer.encode(seq).ids
                else:
                    tokens = self.tokenizer.encode(seq)
                tokenized.append(tokens)
            return tokenized

        # Split sequences into chunks for parallel processing
        chunk_size = len(sequences) // num_workers + 1
        chunks = [
            sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)
        ]

        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            results = list(executor.map(tokenize_batch, chunks))

        # Flatten results
        tokenized_sequences = []
        for chunk_result in results:
            tokenized_sequences.extend(chunk_result)

        return tokenized_sequences

    def group_sequences_by_length(
        self, sequences, names, tokenized_sequences, length_bins=10
    ):
        """Group sequences by similar lengths to minimize padding"""
        seq_lengths = [len(seq) for seq in tokenized_sequences]

        # Create length bins
        min_len, max_len = min(seq_lengths), max(seq_lengths)
        bin_size = (max_len - min_len) // length_bins + 1

        grouped_data = defaultdict(list)

        for i, (seq, name, tokens, length) in enumerate(
            zip(sequences, names, tokenized_sequences, seq_lengths)
        ):
            bin_id = (length - min_len) // bin_size
            grouped_data[bin_id].append((seq, name, tokens, length, i))

        return grouped_data

    def encode_batch_optimized(self, tokenized_seqs, max_batch_length=None):
        """Highly optimized batch encoding"""
        if not tokenized_seqs:
            return []

        # Find max length for minimal padding
        if max_batch_length is None:
            max_len = max(len(seq) for seq in tokenized_seqs)
        else:
            max_len = max_batch_length

        batch_size = len(tokenized_seqs)

        # Pre-allocate tensors for efficiency
        batch_tensor = torch.zeros(
            (batch_size, max_len), dtype=torch.long, device=self.device
        )
        attention_mask = torch.zeros(
            (batch_size, max_len), dtype=torch.long, device=self.device
        )

        # Fill tensors efficiently
        for i, seq in enumerate(tokenized_seqs):
            seq_len = min(len(seq), max_len)
            batch_tensor[i, :seq_len] = torch.tensor(seq[:seq_len], dtype=torch.long)
            attention_mask[i, :seq_len] = 1

        # Forward pass with no gradient computation
        with torch.no_grad():
            logits, hidden_states = self.model(batch_tensor)

        # Efficient mean calculation
        # Multiply by attention mask and divide by sequence lengths
        masked_hidden = hidden_states * attention_mask.unsqueeze(-1).float()
        seq_lengths = attention_mask.sum(dim=1, keepdim=True).float()
        mean_embeddings = masked_hidden.sum(dim=1) / seq_lengths

        return mean_embeddings

    def encode_large_dataset(
        self,
        sequences,
        names,
        batch_size=1000,
        tokenize_workers=4,
        write_buffer_size=10000,
        progress_file=None,
    ):
        """Optimized encoding for very large datasets"""

        print(f"Starting tokenization of {len(sequences)} sequences...")

        # Parallel tokenization
        tokenized_sequences = self.tokenize_sequences_parallel(
            sequences, tokenize_workers
        )

        print("Grouping sequences by length...")
        # Group by length to minimize padding overhead
        grouped_data = self.group_sequences_by_length(
            sequences, names, tokenized_sequences
        )

        results = []
        processed_count = 0

        print("Starting encoding...")

        # Process each length group
        for bin_id in tqdm(sorted(grouped_data.keys()), desc="Processing length bins"):
            bin_data = grouped_data[bin_id]

            # Process this bin in batches
            for i in range(0, len(bin_data), batch_size):
                batch_data = bin_data[i : i + batch_size]
                batch_tokenized = [item[2] for item in batch_data]  # tokens
                batch_names = [item[1] for item in batch_data]  # names
                batch_indices = [item[4] for item in batch_data]  # original indices

                # Encode batch
                embeddings = self.encode_batch_optimized(batch_tokenized)

                # Store results with original indices to maintain order
                for j, (name, embedding, orig_idx) in enumerate(
                    zip(batch_names, embeddings, batch_indices)
                ):
                    results.append((orig_idx, name, embedding.cpu()))

                processed_count += len(batch_data)

                # Memory management
                if processed_count % write_buffer_size == 0:
                    torch.cuda.empty_cache()
                    gc.collect()

                # Save progress periodically
                if progress_file and processed_count % (write_buffer_size * 5) == 0:
                    torch.save(
                        {"processed_count": processed_count, "results": results},
                        progress_file,
                    )

        print(f"Finished encoding {processed_count} sequences")

        # Sort results by original index to maintain input order
        results.sort(key=lambda x: x[0])

        return [(name, embedding) for _, name, embedding in results]

    def write_results_buffered(self, results, output_file, buffer_size=10000):
        """Write results with buffering for better I/O performance"""

        with open(output_file, "w") as f:
            buffer = []

            for name, embedding in tqdm(results, desc="Writing results"):
                out_list = [name] + embedding.tolist()
                line = "\t".join([str(x) for x in out_list]) + "\n"
                buffer.append(line)

                if len(buffer) >= buffer_size:
                    f.writelines(buffer)
                    buffer = []

            # Write remaining buffer
            if buffer:
                f.writelines(buffer)

    def load_model(self, model_cfg, ckpt_path, max_seq_len, nlayer):
        # ... (same as original)
        config = yaml.load(open(model_cfg, "r"), Loader=yaml.FullLoader)
        config["model_config"]["n_layer"] = nlayer
        config["model_config"]["layer"]["l_max"] = max_seq_len + 2
        model = ConvLMHeadModel(**config["model_config"])

        state_dict = torch.load(ckpt_path, map_location="cpu")
        torch.nn.modules.utils.consume_prefix_in_state_dict_if_present(
            state_dict["state_dict"], "model."
        )

        model_state_dict = state_dict["state_dict"]
        for key in list(model_state_dict.keys()):
            if "torchmetrics" in key:
                model_state_dict.pop(key)

        model.load_state_dict(state_dict["state_dict"])

        if config["tokenizer_name"] == "char":
            print("**Using Char-level tokenizer**")
            tokenizer = CharacterTokenizer(
                characters=["A", "C", "G", "T", "N"],
                model_max_length=max_seq_len + 2,
                add_special_tokens=False,
            )
            print(tokenizer._vocab_str_to_int)
        else:
            raise NotImplementedError("You need to provide a custom tokenizer!")

        return model, tokenizer


def load_sequences_memory_efficient(seq_file):
    """Memory efficient sequence loading"""
    sequences = []
    names = []

    with open(seq_file, "r") as f:
        current_seq = ""
        current_name = None

        for line in tqdm(f, desc="Loading sequences"):
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:  # Save previous sequence
                    sequences.append(current_seq)
                    names.append(current_name)
                current_name = line[1:]
                current_seq = ""
            else:
                current_seq += line

        # Don't forget the last sequence
        if current_name is not None:
            sequences.append(current_seq)
            names.append(current_name)

    return sequences, names


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # ... (same arguments as before, plus new ones)
    parser.add_argument("--model_cfg", default="")
    parser.add_argument("--ckpt_path", default="")
    parser.add_argument("--seq_file", default="")
    parser.add_argument("--output_file", default="")
    parser.add_argument("--max_seqlen", default="")
    parser.add_argument("--nlayers", default="")
    parser.add_argument(
        "--batch_size", default="2000", help="Much larger for short sequences"
    )
    parser.add_argument(
        "--tokenize_workers", default="8", help="Parallel tokenization workers"
    )
    parser.add_argument(
        "--write_buffer_size", default="50000", help="Buffer size for writing"
    )
    parser.add_argument("--progress_file", default="", help="File to save progress")
    parser.add_argument("--disable_multi_gpu", action="store_true")

    args = parser.parse_args()

    # Initialize encoder
    use_multi_gpu = not args.disable_multi_gpu
    task = FastHG38Encoder(
        args.model_cfg,
        args.ckpt_path,
        max_seq_len=int(args.max_seqlen),
        nlayer=int(args.nlayers),
        use_multiple_gpus=use_multi_gpu,
    )
    print("Successfully loaded encoder.", flush=True)

    # Load sequences efficiently
    print("Loading sequences...")
    sequences, names = load_sequences_memory_efficient(args.seq_file)

    if len(sequences) != len(names):
        print("Error: Number of sequences and names do not match.", flush=True)
        sys.exit(1)

    print(f"Successfully loaded {len(sequences)} sequences.", flush=True)

    # Process with optimized pipeline
    batch_size = int(args.batch_size)
    tokenize_workers = int(args.tokenize_workers)
    write_buffer_size = int(args.write_buffer_size)

    progress_file = args.progress_file if args.progress_file else None

    # Encode all sequences
    results = task.encode_large_dataset(
        sequences,
        names,
        batch_size=batch_size,
        tokenize_workers=tokenize_workers,
        write_buffer_size=write_buffer_size,
        progress_file=progress_file,
    )

    # Write results efficiently
    print("Writing results to file...")
    task.write_results_buffered(results, args.output_file, write_buffer_size)

    print("Successfully encoded and saved all sequences!", flush=True)
