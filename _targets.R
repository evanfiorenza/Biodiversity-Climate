# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(tidyverse)
library(vegan)
library(betapart)
library(qs)
library(crew)
library(breakaway)

#What happens if I make edits from github on online version of vscode?

#Configure Parallel Controller
controller <- crew::crew_controller_local(
  name = "my_controller",
  workers = 5
)



# Set target options:
tar_option_set(
  packages = c("tidyverse","vegan","betapart","glmmTMB","breakaway","crew","qs"),
  memory = "transient",
  garbage_collection = TRUE,
  format="qs",
  error="null",
  controller=controller,
  storage = "worker",
  retrieval = "worker"
  
  
  
   # Packages that your targets need for their tasks.
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source('R/functions.R')
# tar_source("other_functions.R") # Source other scripts as needed.
# Replace the target list below with your own:
list(
  tar_target(
    name=BioTime.data,
    command = get_data_csv("Data/Biotime_Renamed_Gridded_12_Filtered.csv")
   ),
  tar_target(
    name=BioTime.data2,
    command = group_data_csv(BioTime.data)
  ),
  mapped<-tar_map(
    values=get_data_csv("Data/Values_AssemblageID.csv"),
    tar_target(
    name=spp.matrix,
    command = clean_and_pivot(BioTime.data2,"Species","ABUNDANCE",0,assemblage.ID),
   
   ),
 #  tar_target(
  #   name=alpha.table,
  #  command=get_alpha_div_metrics(spp.matrix,10)
  # ),
  # tar_target(
  #   name=alpha.trends,
  #   command=get_alpha_div_trends(alpha.table)
  # ),
  tar_target(
    name=decomposed,
    command = get_decomposed_beta_pairwise(spp.matrix,index.family = "jaccard",10),
    
  ),
  tar_target(
    name=extracted.beta,
    command = extract_decomposed_beta_start_to_end(decomposed)
  )
),
tar_combine(
  name = combined.beta,
  mapped[['extracted.beta']]
)#,
#tar_combine(
#  name = combined.alpha,
#  mapped[['alpha.trends']]
#)
)






