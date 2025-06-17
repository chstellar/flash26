library(tidyverse)

dir = "/oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline/results/"
all_files <- list.files(path = dir, pattern = "_blast_annotated_plots.pdf|_heatmaps.pdf|confusion_matrices.pdf", 
                        recursive = T, full.names = T)

dir = "/oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline/results//"

out_dir <- "/oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline/results/summary/250610/"

dir.create(out_dir, recursive = T)

my_files <- data.frame(file=all_files)

my_files <- my_files %>% mutate(new_file=all_files %>% 
                                  str_remove(dir) %>% 
                                  str_remove("(?<=ized\\/).+(?=adelie)") %>%
                                  str_replace_all("/", "_")) %>% 
  mutate(new_dir = str_extract(new_file, "(.+)_filter", group=1))

my_files <- my_files %>% filter(str_detect(new_file, "filter1")) %>% 
  filter(str_detect(file, "top20000")) %>%
  filter(str_detect(new_file, "ohe|hyena_nor")) %>% 
  filter(str_detect(new_file, "shiftDist-lev"))

my_files <- my_files %>% filter(!str_detect(file, "train_test_split|/old/"))

my_files <- my_files %>% filter(!str_detect(file, "summary")) %>% filter(!str_detect(file, "withinOrder"))

new_dirs <- unique(my_files$new_dir)

map(file.path(out_dir, new_dirs), dir.create, recursive=TRUE)

my_files <- my_files %>% mutate(out_file = file.path(out_dir, new_dir, new_file))

walk2(my_files$file, my_files$out_file, \(x,y) system(paste("cp", x, y)))

### pneumo 

my_files <- data.frame(file=all_files)

my_files <- my_files %>% mutate(new_file=all_files %>% 
                                  str_remove(dir) %>% 
                                  str_remove("(?<=ized\\/).+(?=adelie)") %>%
                                  str_replace_all("/", "_")) %>% 
  mutate(new_dir = str_extract(new_file, "(.+)_filter", group=1))

my_files <- my_files %>% filter(str_detect(new_file, "filter1")) %>% 
  filter(str_detect(file, "top20000")) %>%
  filter(str_detect(new_file, "pneumo")) %>%
  filter(str_detect(new_file, "hyena_unn")) %>% 
  filter(str_detect(new_file, "shiftDist-lev|masked"))

my_files <- my_files %>% filter(!str_detect(file, "train_test_split|/old/|summary"))

my_files <- my_files %>% mutate(new_dir = paste0(new_dir, "_hyena_unnormalized")) 
my_files <- my_files %>% filter(!str_detect(file, "masked-nucleotide"))

new_dirs <- unique(my_files$new_dir)

map(file.path(out_dir, new_dirs), dir.create, recursive=TRUE)

my_files <- my_files %>% mutate(out_file = file.path(out_dir, new_dir, new_file))
walk2(my_files$file, my_files$out_file, \(x,y) system(paste("cp", x, y)))

### y1000 - masked aa 


my_files <- data.frame(file=all_files)

my_files <- my_files %>% mutate(new_file=all_files %>% 
                                  str_remove(dir) %>% 
                                  str_remove("(?<=ized\\/).+(?=adelie)") %>%
                                  str_replace_all("/", "_")) %>% 
  mutate(new_dir = str_extract(new_file, "(.+)_filter", group=1))

my_files <- my_files %>% filter(str_detect(new_file, "filter1")) %>% 
  filter(str_detect(file, "top20000")) %>%
  filter(str_detect(new_file, "y1000-genomes-resistance|y1000-genomes-data")) %>%
  filter(str_detect(new_file, "hyena_nor")) %>% 
  filter(str_detect(new_file, "masked"))

my_files <- my_files %>% filter(!str_detect(file, "train_test_split|/old/|summary"))

my_files <- my_files %>% mutate(new_dir = paste0(new_dir, "_masked-aa-clustered")) 
my_files <- my_files %>% filter(!str_detect(file, "masked-nucleotide"))

new_dirs <- unique(my_files$new_dir)

map(file.path(out_dir, new_dirs), dir.create, recursive=TRUE)

my_files <- my_files %>% mutate(out_file = file.path(out_dir, new_dir, new_file))
walk2(my_files$file, my_files$out_file, \(x,y) system(paste("cp", x, y)))

