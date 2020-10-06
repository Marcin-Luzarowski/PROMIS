require(tidyverse)
require(pastecs)
require(formattable)
require(pheatmap)
require(RColorBrewer)
require(gplots)

#Before start - load functions, which are present in:
#1) load_data.R
#2) plot_data.R
#3) data_processing.R
#4) check training data
#NOTE: Description of all the function, required input and what function return can be find in separate source files.

# set working directory
wd <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example"
setwd(wd)

#Provide pathway to your protein input file
protein_input <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Protein_example.txt"
metabolite_input <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Metabolite_example.txt"

###########################################################################################
#Import protein dataset
#it requires to provide number of replicates in protein dataset
###########################################################################################
#Protein data
nr_replicas <- 1
protein_df <- import_data(protein_input, nr_replicas = nr_replicas)

#Metabolite data
nr_replicas <- 3
metabolite_df <- import_data(metabolite_input, nr_replicas = nr_replicas)

###########################################################################################
#Remove single appearances and short peaks#
#it requires to specify number of collected fractions and number of replicas
#input is either protein or metabolite dataframe
###########################################################################################

#Protein data
nr_replicas <- 1
nr_fractions <- 38
protein_long <- rmv_short_peaks(protein_df, nr_replicas = nr_replicas, nr_fractions = nr_fractions)

#Metabolite data
nr_replicas <- 3
nr_fractions <- 38
metabolite_long <- rmv_short_peaks(metabolite_df, nr_replicas = nr_replicas, nr_fractions = nr_fractions)


###########################################################################################
#Select molecules having reproducible elution profiles#
#the input is preferably dataframe with removed short peaks
#one has to define required number of replicates having reproducible elution profile
#one has to define minimal similarity level between the profiles
###########################################################################################

#Protein data
req_repl <- 1
repr_profile_threshold <- 0.0

protein_reproducible <- repr_profile(protein_long, req_repl = req_repl, repr_profile_threshold = repr_profile_threshold)

#Metabolite data
req_repl <- 2
repr_profile_threshold <- 0.7

met_reproducible <- repr_profile(metabolite_long, req_repl = req_repl, repr_profile_threshold = repr_profile_threshold)

# Calculate single profile using data from x replicas
# define singe_profile = "sum" or "mean" or "median"
# define nr_replicas

met_single_profile <- calc_single_profile(x = met_reproducible, single_profile = "median", nr_replicas = 3)
prot_single_profile <- calc_single_profile(x = protein_reproducible, single_profile = "median", nr_replicas = 1)

#' Plot elution profile of a molecule
#' @param x Data frame containing single elution profiles
#' @param met_plot_color Color of a plot
#' @param do_plot If plot should be made
#' @param normalize If data should be maximum transformed
#' @param tmp_dir Determine temporary directory to which the peaks will be saved
#' @return A data frame without elution profiles, which peaks span less than c fractions
#' @export

# Plot single profiles of metabolites
met_for_deconv <- plot_single_profiles(x = met_single_profile,
                                       met_plot_color = "blue",
                                       do_plot = TRUE,
                                       normalize = TRUE,
                                       tmp_dir = "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Plots/Metabolites",
                                       fraction_names = paste("log_", c(1:nr_fractions), sep = ""))

# Plot single profiles of proteins

prot_for_deconv <- plot_single_profiles(x = prot_single_profile,
                                        met_plot_color = "red",
                                        do_plot = TRUE,
                                        normalize = TRUE,
                                        tmp_dir = "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Plots/Proteins",
                                        fraction_names = paste("log_", c(1:nr_fractions), sep = ""))

# Deconvolute elution profiles of metabolites
met_deconvoluted_data <- deconvolution(x = met_for_deconv,
                                       var_min_peak = 0.2,
                                       var_limit = 0.15,
                                       var_limit_down = 0.75,
                                       var_limit_up = 0.75,
                                       var_num_col = 38,
                                       var_size_range_name = "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Column calibration.txt",
                                       var_exp_name = paste("log_deconvoluted_metabolites"),
                                       var_work_table = 3,
                                       var_plot_decon_peaks = 2,
                                       dec_data_color = "grey",
                                       plot_single_plot = FALSE,
                                       tmp_dir = "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Plots/Deconvoluted/",
                                       fraction_names = paste(c(1:nr_fractions), sep = ""))

# Deconvolute elution profiles of proteins
prot_deconvoluted_data <- deconvolution(x = prot_for_deconv,
                                        var_min_peak = 0.2,
                                        var_limit = 0.15,
                                        var_limit_down = 0.75,
                                        var_limit_up = 0.75,
                                        var_num_col = 38,
                                        var_size_range_name = "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Column calibration.txt",
                                        var_exp_name = paste("log_deconvoluted_proteins"),
                                        var_work_table = 3,
                                        var_plot_decon_peaks = 2,
                                        dec_data_color = "grey",
                                        plot_single_plot = FALSE,
                                        tmp_dir = "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Plots/Deconvoluted/",
                                        fraction_names = paste(c(1:nr_fractions), sep = ""))

# Perform adduct detection
df_adduct_detection <- read.delim("Adduct_detection_example.txt")
df_adducts <- read.delim("Adduct_pos.txt")
test_adduct <- adduct_detection(x = df_adduct_detection,
                                y = df_adducts,
                                mass_deviation = 0.015,
                                RT_deviation = 0.005)

# Prepare correlation table
correlation_matrix <- calc_correlation(x = prot_deconvoluted_data,
                                       y = met_deconvoluted_data)

# Find unique features
FUF_df <- read.delim("FUF_PCF.txt")
FUF_info <- read.delim("FUF_info.txt")

unique_features <- find_unique_features(FUF_df = FUF_df, FUF_info = FUF_info, nr_fractions = 38,
                                        nr_replicas = 3, start_column = 1, end_column = 114)

# Find proteins eluting in a complex
deconv_profiles <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Prot_dec.txt"
prot_theor_mass <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Yeast_protein.txt"
column_calibration <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Mass_fractions.txt"
oligomeric_state_ratio <- 1.5

proteins_in_complex <- find_proteins_in_complex(deconv_profiles = deconv_profiles,
                                                prot_theor_mass = prot_theor_mass,
                                                column_calibration = column_calibration,
                                                oligomeric_state_ratio = oligomeric_state_ratio)

# Find elution profiles of known protein-protein complexes
deconv_profiles <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Prot_dec.txt"
prot_theor_mass <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Yeast_protein.txt"
column_calibration <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Mass_fractions.txt"
oligomeric_state_ratio <- 1.5
def_heatmap_name <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/PPI_Known_test.jpeg"
ref_database <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/4932.protein.links.full.v10.5.txt"

known_PPI_heatmap <- find_known_PPI(deconv_profiles = deconv_profiles,
                                    prot_theor_mass = prot_theor_mass,
                                    column_calibration = column_calibration,
                                    oligomeric_state_ratio = oligomeric_state_ratio,
                                    ref_database = ref_database,
                                    def_heatmap_name = def_heatmap_name)

# Find known protein-metabolite complexes and prepare ROC curve
ref_library <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/4932.protein_chemical.links.detailed.v5.0.tsv"
prot_theor_mass <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Yeast_protein.txt"
prot_db <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Filt_prot_merged.txt"
annotated_met <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/MetabolitesInPCF_no_red.txt"
met_chem_id <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/MetaboliteLibrary_ChemForm_CID.txt"
cor_matrix <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Corr_prot_vs_sm.txt"

known_PMI <- find_known_PMI(ref_library = ref_library,
                            prot_theor_mass = prot_theor_mass,
                            prot_db = prot_db,
                            annotated_met = annotated_met,
                            met_chem_id = met_chem_id,
                            cor_matrix = cor_matrix)

# Find predicted protein-metabolite complexes
ref_library <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/4932.protein_chemical.links.detailed.v5.0.tsv"
prot_theor_mass <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Yeast_protein.txt"
prot_db <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Filt_prot_merged.txt"
annotated_met <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/MetabolitesInPCF_no_red.txt"
met_chem_id <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/MetaboliteLibrary_ChemForm_CID.txt"
cor_matrix <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Corr_prot_vs_sm.txt"

predicted_PMI <- find_predicted_PMI(ref_library = ref_library,
                                    prot_theor_mass = prot_theor_mass,
                                    prot_db = prot_db,
                                    annotated_met = annotated_met,
                                    met_chem_id = met_chem_id,
                                    cor_matrix = cor_matrix)

# Determine accumulation pattern of metabolites
metabolites_df <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/peaks_pp.txt"
annotation_df <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/annotation_DP.txt"
s_order <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/DP_accumulation_order.txt"
tmp_file_name <- "Dipeptide_accumulation.jpeg"

dipeptides_example <- metabolites_accumulation(metabolites_df = metabolites_df,
                                               annotation_df = annotation_df,
                                               s_order = s_order,
                                               tmp_file_name = tmp_file_name)

# Determine 13C fraction size
enrichment_df <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/enrichment_dummy.txt"
s_order <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/enrichment_order.txt"
metabolite_poolsize <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/met_poolsize_dummy.txt"
lab_annotation <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/13C_annotation.txt"
tmp_barplot_name <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/Barplots.pdf"
tmp_dir_plot <- "H:/PhD/Projects/4. Yeast/!!!Yeast_PROMIS/!R Project/!2020_Submission to GitHub/Data_example/13C_Plots/"

signif_changes_13C <- analyze_13C_data(enrichment_df = enrichment_df,
                                       s_order = s_order,
                                       metabolite_poolsize = metabolite_poolsize,
                                       lab_annotation = lab_annotation,
                                       tmp_barplot_name = tmp_barplot_name,
                                       tmp_dir_plot = tmp_dir_plot)
