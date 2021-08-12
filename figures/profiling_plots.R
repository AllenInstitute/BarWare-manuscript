library(tidyverse)
library(janitor)
library(gridExtra)
library(cowplot)

setwd("G:\\Shared drives\\Imm - All\\Manuscripts\\2021 Cell Hashing Pipeline\\ADT_Count_Profiling")

results <- read_csv("ADT_profiling_results.csv")

# clean data
results <- janitor::clean_names(results)

results <- results %>%
  # Drop non-hashed sample
  filter(well != "NON-1-HT5") %>%
  # Add N cells per well as a column
  mutate(well_n = as.integer(sub("Pool-","",well))) %>%
  # cast percent CPUs as int
  mutate(percent_cpu = as.integer(percent_cpu)) %>%
  # convert wall time to s
  mutate(elapsed_time_s = as.integer(elapsed_time_h_mm_ss))


# original version ----------------------------

# memory
mem <- results %>%
  ggplot(aes(x = well, y = max_resident_set_size_kb, color = tool)) +
  geom_point() +
  theme_bw() +
  labs(title = "Maximum Memory in KB (resident set size)",
       x = "sample well",
       y = "maximum resident set size (KB)")

# user time
u_time <- results %>%
  ggplot(aes(x = well, y = user_time_s, color = tool)) +
  geom_point() +
  theme_bw() +
  labs(title = "User Time in Seconds",
       x = "sample well",
       y = "user time (seconds)")

# elapsed (wall clock) time
wall_time <- results %>%
  ggplot(aes(x = well, y = elapsed_time_s, color = tool)) +
  geom_point() +
  theme_bw() +
  labs(title = "Elapsed (Wall Clock) Time",
       x = "sample well",
       y = "elapsed time (hr:min:sec)")

# percent CPU
cpu <- results %>%
  ggplot(aes(x = well, y = percent_cpu, color = tool)) +
  geom_point() +
  theme_bw() +
  labs(title = "Average Percent CPU Used",
       x = "sample well",
       y = "Percent CPU")

# output size
log_coutput_size <- results %>%
  ggplot(aes(x = well, y = log10(output_size_bytes), color = tool)) +
  geom_point() +
  theme_bw() +
  labs(title = "Output Size in Bytes (log scale)",
       x = "N Cells loaded",
       y = "log10 of output size in bytaes")


## Lucas' version ---------------------------

## memory
mem2 <- ggplot() +
  geom_line(data = results,
            aes(x = well_n,
                y = max_resident_set_size_kb / 1e6,
                group = tool,
                color = tool),
            size = 1) +
  geom_point(data = results,
             aes(x = well_n,
                 y = max_resident_set_size_kb / 1e6,
                 group = tool,
                 color = tool),
             size = 2) +
  theme_bw() +
  scale_x_continuous(breaks = c(16, 24, 32, 48, 64, 80)) +
  scale_y_continuous(limits = c(0, 24), breaks = seq(0,24,4)) +
  labs(title = "Maximum Memory in GB (resident set size)",
       x = "N Cells loaded",
       y = "maximum resident set size (GB)") +
  theme(legend.position = "none")

## elapsed (wall clock) time
wall_time2 <- ggplot() +
  geom_line(data = results,
            aes(x = well_n,
                y = elapsed_time_s / 60,
                group = tool,
                color = tool),
            size = 1) +
  geom_point(data = results,
             aes(x = well_n,
                 y = elapsed_time_s / 60,
                 group = tool,
                 color = tool),
             size = 2) +
  theme_bw() +
  scale_x_continuous(breaks = c(16, 24, 32, 48, 64, 80)) +
  scale_y_log10(limits = c(1,500), breaks = c(1, 5, 15, 60, 120, 240, 480)) +
  labs(title = "Elapsed (Wall Clock) Time",
       x = "N Cells loaded",
       y = "elapsed time (minutes)")

## user time
u_time2 <- ggplot() +
  geom_line(data = results,
            aes(x = well_n,
                y = user_time_s / 60,
                group = tool,
                color = tool),
            size = 1) +
  geom_point(data = results,
             aes(x = well_n,
                 y = user_time_s / 60,
                 group = tool,
                 color = tool),
             size = 2) +
  theme_bw() +
  scale_x_continuous(breaks = c(16, 24, 32, 48, 64, 80)) +
  scale_y_log10(limits = c(1,500), breaks = c(1, 5, 15, 60, 120, 240, 480)) +
  labs(title = "User Time in Minutes",
       x = "sample well",
       y = "user time (minutes)")

## percent CPU
cpu2 <- ggplot() +
  geom_line(data = results,
            aes(x = well_n,
                y = percent_cpu,
                group = tool,
                color = tool),
            size = 1) +
  geom_point(data = results,
             aes(x = well_n,
                 y = percent_cpu,
                 group = tool,
                 color = tool),
             size = 2) +
  theme_bw() +
  scale_x_continuous(breaks = c(16, 24, 32, 48, 64, 80)) +
  scale_y_continuous(limits = c(0,600), breaks = seq(0, 600, 100)) +
  labs(title = "Average Percent CPU Load",
       x = "N Cells loaded",
       y = "percent CPU load") +
  theme(legend.position = "none")


# merge plots ---------------------------
gridExtra::grid.arrange(u_time, mem, cpu, log_coutput_size)

plot_grid(mem2, cpu2, wall_time2, 
          rel_widths = c(1,1,1.6),
          ncol = 3)

