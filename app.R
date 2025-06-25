options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("shiny")
library(shiny)
library(shinyjs)
library(rsconnect)


# rsconnect::writeManifest( # Commented out as this is for deployment, not app logic
#   appDir = ".",
#   appPrimaryDoc = "app.R",
#   appFiles = NULL
# )

# Define the genetic code (mRNA codons to amino acids)
genetic_code <- list(
  "UUU" = "Phe", "UUC" = "Phe", "UUA" = "Leu", "UUG" = "Leu",
  "CUU" = "Leu", "CUC" = "Leu", "CUA" = "Leu", "CUG" = "Leu",
  "AUU" = "Ile", "AUC" = "Ile", "AUA" = "Ile", "AUG" = "Met", # AUG is also Start
  "GUU" = "Val", "GUC" = "Val", "GUA" = "Val", "GUG" = "Val",
  "UCU" = "Ser", "UCC" = "Ser", "UCA" = "Ser", "UCG" = "Ser",
  "CCU" = "Pro", "CCC" = "Pro", "CCA" = "Pro", "CCG" = "Pro",
  "ACU" = "Thr", "ACC" = "Thr", "ACA" = "Thr", "ACG" = "Thr",
  "GCU" = "Ala", "GCC" = "Ala", "GCA" = "Ala", "GCG" = "Ala",
  "UAU" = "Tyr", "UAC" = "Tyr", "UAA" = "STOP", "UAG" = "STOP",
  "CAU" = "His", "CAC" = "His", "CAA" = "Gln", "CAG" = "Gln",
  "AAU" = "Asn", "AAC" = "Asn", "AAA" = "Lys", "AAG" = "Lys",
  "GAU" = "Asp", "GAC" = "Asp", "GAA" = "Glu", "GAG" = "Glu",
  "UGU" = "Cys", "UGC" = "Cys", "UGA" = "STOP", "UGG" = "Trp",
  "CGU" = "Arg", "CGC" = "Arg", "CGA" = "Arg", "CGG" = "Arg",
  "AGU" = "Ser", "AGC" = "Ser", "AGA" = "Arg", "AGG" = "Arg",
  "GGU" = "Gly", "GGC" = "Gly", "GGA" = "Gly", "GGG" = "Gly"
)

# List of valid 3-letter amino acid codes
valid_amino_acid_codes <- toupper(unlist(unname(genetic_code)))
valid_amino_acid_codes <- unique(valid_amino_acid_codes[valid_amino_acid_codes != "STOP"])
valid_amino_acid_codes <- c(valid_amino_acid_codes, "STOP")

# Define amino acid categories and their colors based on the provided image
amino_acid_properties <- list(
  # Aliphatic (Red)
  "ALA" = list(category = "Aliphatic", color = "#FF0000"), # Alanine - Red
  "GLY" = list(category = "Aliphatic", color = "#FF0000"), # Glycine - Red
  "ILE" = list(category = "Aliphatic", color = "#FF0000"), # Isoleucine - Red
  "LEU" = list(category = "Aliphatic", color = "#FF0000"), # Leucine - Red
  "PRO" = list(category = "Aliphatic", color = "#FF0000"), # Proline - Red
  "VAL" = list(category = "Aliphatic", color = "#FF0000"), # Valine - Red
  
  # Aromatic (Green)
  "PHE" = list(category = "Aromatic", color = "#99CC66"), # Phenylalanine - Light Green
  "TRP" = list(category = "Aromatic", color = "#99CC66"), # Tryptophan - Light Green
  "TYR" = list(category = "Aromatic", color = "#99CC66"), # Tyrosine - Light Green
  
  # Acidic (Deep Orange)
  "ASP" = list(category = "Acidic", color = "#FF6F00"), # Aspartic Acid - Deep Orange
  "GLU" = list(category = "Acidic", color = "#FF6F00"), # Glutamic Acid - Deep Orange
  
  # Basic (Punchy Cyan Blue)
  "ARG" = list(category = "Basic", color = "#00B7EB"), # Arginine - Punchy Cyan Blue
  "HIS" = list(category = "Basic", color = "#00B7EB"), # Histidine - Punchy Cyan Blue
  "LYS" = list(category = "Basic", color = "#00B7EB"), # Lysine - Punchy Cyan Blue
  
  # Hydroxylic (Pink)
  "SER" = list(category = "Hydroxylic", color = "#FFB6D5"), # Serine - Pink
  "THR" = list(category = "Hydroxylic", color = "#FFB6D5"), # Threonine - Pink
  
  # Sulfur-containing (Vivid Yellow/Gold)
  "CYS" = list(category = "Sulfur-containing", color = "#FFD700"), # Cysteine - Vivid Yellow/Gold
  "MET" = list(category = "Sulfur-containing", color = "#FFD700"), # Methionine - Vivid Yellow/Gold
  
  # Amidic (Dark Blue)
  "ASN" = list(category = "Amidic", color = "#003f88"), # Asparagine - Dark Blue
  "GLN" = list(category = "Amidic", color = "#003f88"), # Glutamine - Dark Blue
  
  # Special case for STOP codon
  "STOP" = list(category = "Stop Codon", color = "#696969") # Dim Gray for stop
)

# Function to validate DNA sequence
validate_dna <- function(dna) {
  dna <- toupper(gsub(" ", "", dna))
  if (!grepl("^[ATCG]*$", dna)) {
    return("Error: DNA sequence must contain only A, T, C, G")
  }
  return(dna)
}

# Color mapping for amino acids (3-letter codes) based on codon chart colors
# CSS Class Suggestions for future use:
# .aa-basic { color: #00B7EB; }
# .aa-amidic { color: #003f88; }
# .aa-acidic { color: #FF6F00; }
# .aa-sulfur { color: #FFD700; }
amino_acid_properties <- list(
  # Aliphatic (Red)
  "ALA" = list(category = "Aliphatic", color = "#FF0000"),
  "GLY" = list(category = "Aliphatic", color = "#FF0000"),
  "ILE" = list(category = "Aliphatic", color = "#FF0000"),
  "LEU" = list(category = "Aliphatic", color = "#FF0000"),
  "PRO" = list(category = "Aliphatic", color = "#FF0000"),
  "VAL" = list(category = "Aliphatic", color = "#FF0000"),

  # Aromatic (Green)
  "PHE" = list(category = "Aromatic", color = "#99CC66"),
  "TRP" = list(category = "Aromatic", color = "#99CC66"),
  "TYR" = list(category = "Aromatic", color = "#99CC66"),

  # Acidic (Deep Orange)
  "ASP" = list(category = "Acidic", color = "#FF6F00"),
  "GLU" = list(category = "Acidic", color = "#FF6F00"),

  # Basic (Punchy Cyan Blue)
  "ARG" = list(category = "Basic", color = "#00B7EB"),
  "HIS" = list(category = "Basic", color = "#00B7EB"),
  "LYS" = list(category = "Basic", color = "#00B7EB"),

  # Hydroxylic (Pink)
  "SER" = list(category = "Hydroxylic", color = "#FFB6D5"),
  "THR" = list(category = "Hydroxylic", color = "#FFB6D5"),

  # Sulfur-containing (Vivid Yellow/Gold)
  "CYS" = list(category = "Sulfur-containing", color = "#FFD700"),
  "MET" = list(category = "Sulfur-containing", color = "#FFD700"),

  # Amidic (Dark Blue)
  "ASN" = list(category = "Amidic", color = "#003f88"),
  "GLN" = list(category = "Amidic", color = "#003f88"),

  # Special case for STOP codon
  "STOP" = list(category = "Stop Codon", color = "#696969")
)

# Function to format amino acid sequence with colors based on properties
format_amino_acids_with_color <- function(amino_acid_string) {
  if (is.null(amino_acid_string) || nchar(amino_acid_string) == 0) {
    return(NULL) # Return NULL for empty, so it doesn't render empty spans
  }
  
  amino_acid_string <- toupper(gsub(" ", "", amino_acid_string))
  amino_acids <- strsplit(amino_acid_string, "(?<=\\G...)", perl = TRUE)[[1]]
  
  formatted_html_elements <- lapply(amino_acids, function(aa_code) {
    props <- amino_acid_properties[[aa_code]]
    if (is.null(props)) {
      return(tags$span(class = "colored-amino-acid", style = "color: #333333; font-weight: bold;", aa_code)) # Add bold
    } else {
      return(tags$span(class = "colored-amino-acid", style = paste0("color: ", props$color, "; font-weight: bold;"), aa_code)) # Add bold
    }
  })
  
  # Simply combine the individual span elements. CSS will handle spacing.
  # Use HTML() to ensure the list of tags is rendered as HTML, rather than R's default list printing.
  # You might need to flatten the list if it's deeply nested, but `tagList` usually handles this.
  tagList(formatted_html_elements)
}

# Function to transcribe DNA to RNA
transcribe_dna <- function(dna) {
  dna <- toupper(gsub(" ", "", dna))
  rna <- chartr("ATCG", "UAGC", dna)
  return(rna)
}

# Function to validate RNA sequence (no longer user input, but still needs format check)
validate_rna_format <- function(rna_seq) {
  rna_seq <- toupper(gsub(" ", "", rna_seq))
  if (!grepl("^[AUCG]*$", rna_seq)) {
    return("Error: RNA sequence contains invalid characters. Only A, U, C, G are allowed.")
  }
  if (nchar(rna_seq) %% 3 != 0) {
    return("Error: RNA sequence length must be a multiple of 3 for translation.")
  }
  return(rna_seq)
}

# Function to translate RNA to Amino Acids
translate_rna <- function(rna_seq) {
  rna_seq <- toupper(gsub(" ", "", rna_seq))
  if (nchar(rna_seq) %% 3 != 0) {
    return("Error: RNA sequence length must be a multiple of 3 for translation.")
  }
  codons <- strsplit(rna_seq, "(?<=\\G...)", perl = TRUE)[[1]]
  amino_acids <- sapply(codons, function(codon) {
    aa <- genetic_code[[codon]]
    if (is.null(aa)) {
      return("Invalid Codon") 
    }
    return(aa)
  })
  return(amino_acids)
}

# Function to validate Amino Acid sequence entered by user
validate_amino_acids_input <- function(expected_rna_seq, user_amino_acid_input) {
  user_amino_acid_input <- toupper(gsub(" ", "", user_amino_acid_input))
  
  # Break user input into 3-letter pieces
  input_aas_chunks <- strsplit(user_amino_acid_input, "(?<=\\G...)", perl = TRUE)[[1]]
  
  # Check that 3-letter pieces are all valid amino acid codes
  if (any(!input_aas_chunks %in% valid_amino_acid_codes)) {
    invalid_codes <- input_aas_chunks[!input_aas_chunks %in% valid_amino_acid_codes]
    # Add spaces between codes for clearer error message
    invalid_codes_str <- paste(invalid_codes, collapse = ", ")
    return(paste("Error: Invalid amino acid code(s) entered:", invalid_codes_str, ". Please use standard 3-letter codes (e.g., MET, TRP, STOP)."))
  }
  
  expected_amino_acids_array <- translate_rna(expected_rna_seq)
  
  # Handle potential errors from translate_rna
  if (grepl("^Error:", expected_amino_acids_array[1])) {
    return(expected_amino_acids_array[1])
  }
  
  expected_amino_acids_string <- paste(expected_amino_acids_array, collapse = "")
  expected_amino_acids_string <- toupper(expected_amino_acids_string)
  
  # Check that the length of the amino acid sequence is correct
  if (nchar(user_amino_acid_input) != nchar(expected_amino_acids_string)) {
    return("Error: Amino acid sequence length is incorrect.")
  }
  
  # Check that the amino acid sequence is a match based the amino acid codon chart
  if (user_amino_acid_input != expected_amino_acids_string) {
    return("Error: Amino acid sequence does not match the translated RNA sequence.")
  }
  
  return(user_amino_acid_input)
}

# Function to format DNA/RNA sequence for display (just adds spaces)
format_sequence <- function(seq_string, chunk_size = 3) {
  seq_string <- as.character(seq_string) 
  # Check if the sequence is empty to avoid errors with substring
  if (nchar(seq_string) == 0) {
    return("")
  }
  # Add spaces between chunks but limit the total length
  formatted <- paste(substring(seq_string, seq(1, nchar(seq_string), chunk_size), seq(chunk_size, nchar(seq_string), chunk_size)), collapse = " ")
  # If the formatted string is too long, truncate it
  if (nchar(formatted) > 100) {
    formatted <- substr(formatted, 1, 100)
    formatted <- paste0(formatted, "...")
  }
  return(formatted)
}

# Function to generate random DNA sequence
generate_dna <- function(length = 15) {
  # Always use 15 nucleotides for consistency
  length <- 15
  
  # Function to check if a DNA sequence would create stop codons when transcribed
  has_stop_codons <- function(dna_seq) {
    # Check for stop codons: UAA, UAG, UGA
    # These correspond to DNA sequences: TAA, TAG, TGA
    stop_patterns <- c("TAA", "TAG", "TGA")
    
    for (pattern in stop_patterns) {
      if (grepl(pattern, dna_seq)) {
        return(TRUE)
      }
    }
    return(FALSE)
  }
  
  # Generate DNA sequence and check for stop codons
  max_attempts <- 500  # Increased attempts for better safety
  for (attempt in 1:max_attempts) {
    # Use only safe bases to avoid stop codons
    safe_bases <- c("A", "C", "G")  # Avoid T completely to prevent stop codons
    dna_seq <- paste(sample(safe_bases, length, replace = TRUE), collapse = "")
    
    # Double-check by transcribing and verifying no stop codons in RNA
    rna_seq <- chartr("ATCG", "UAGC", dna_seq)
    if (!grepl("UAA|UAG|UGA", rna_seq)) {
      return(dna_seq)
    }
  }
  
  # If we still can't find a safe sequence, generate a completely safe one
  safe_bases <- c("A", "C", "G")
  return(paste(sample(safe_bases, length, replace = TRUE), collapse = ""))
}

# UI Definition
ui <- fluidPage(
  useShinyjs(),
  titlePanel("🌟 Central Dogma Explorer"),
  tags$style(HTML("
    body {
      background-color: #f0f4f8;
      font-family: 'Inter', sans-serif;
      color: #333;
      overflow-x: hidden;
      margin: 0;
      padding: 0;
    }
    .title-panel {
      text-align: center;
      color: #2c3e50;
      font-size: 28px;
      margin-bottom: 20px;
    }
    .intro-text {
      text-align: center;
      font-size: 16px;
      color: #666;
      margin-bottom: 30px;
      background-color: #e6f0fa;
      padding: 15px;
      border-radius: 8px;
    }
    .input-section {
      background-color: #ffffff;
      padding: 20px;
      border-radius: 10px;
      box-shadow: 0 4px 10px rgba(0,0,0,0.1);
      margin-bottom: 20px;
      width: 100%;
      box-sizing: border-box;
    }
    .input-field {
      margin-bottom: 15px;
    }
    .input-field label {
      font-weight: bold;
      font-size: 18px;
      color: #2c3e50;
      display: block;
      margin-bottom: 5px;
    }
    .input-field input {
      font-family: 'Courier New', monospace !important;
      font-size: 22px !important;
      padding: 12px !important;
      border-radius: 8px !important;
      border: 2px solid #ccc !important;
      width: 100% !important;
      height: 60px !important;
      box-sizing: border-box !important;
      overflow-x: hidden !important;
      white-space: normal !important;
      word-wrap: break-word !important;
      word-break: break-all !important;
      min-width: 0 !important;
      max-width: none !important;
    }
    .input-field input:focus {
      border-color: #007bff !important;
      outline: none !important;
    }
    .valid {
      border-color: #28a745 !important;
    }
    .invalid {
      border-color: #dc3545 !important;
    }
    .button-section {
      text-align: center;
      margin: 20px 0;
    }
    .btn-action {
      font-size: 18px;
      padding: 12px 24px;
      border-radius: 5px;
      margin: 0 10px 10px 0;
      transition: background-color 0.3s ease-in-out, transform 0.2s ease-in-out;
      box-shadow: 0 4px 6px rgba(0,0,0,0.1);
      border: none;
    }
    .btn-process {
      background-color: #007bff;
      color: white;
    }
    .btn-process:hover {
      background-color: #0056b3;
    }
    .btn-refresh {
      background-color: #28a745;
      color: white;
    }
    .btn-refresh:hover {
      background-color: #218838;
    }
    .results-section {
      background-color: #ffffff;
      padding: 20px;
      border-radius: 10px;
      box-shadow: 0 4px 10px rgba(0,0,0,0.1);
      width: 100%;
      box-sizing: border-box;
    }
    .sequence-label {
      font-weight: bold;
      font-size: 20px;
      color: #2c3e50;
      margin-top: 15px;
      margin-bottom: 8px;
    }
    .sequence-text {
      font-family: 'Courier New', monospace;
      font-size: 20px;
      line-height: 1.3;
      background-color: #f8f9fa;
      padding: 15px;
      border-radius: 5px;
      border: 1px solid #ddd;
      min-height: 60px;
      word-wrap: break-word;
      word-break: break-all;
      white-space: normal;
      width: 100%;
      box-sizing: border-box;
      overflow: hidden;
    }
    .sequence-text span.colored-amino-acid {
      margin-right: 10px;
      margin-bottom: 6px;
      font-family: 'Courier New', monospace;
      font-size: 20px;
      font-weight: bold;
      display: inline-block;
      min-width: 50px;
      text-align: center;
      padding: 4px 8px;
      border-radius: 4px;
      background-color: #ffffff;
    }
    .error-message {
      color: #dc3545;
      font-size: 18px;
      margin-top: 10px;
      font-weight: bold;
    }
    .success-message {
      color: #28a745;
      font-size: 18px;
      margin-top: 10px;
      font-weight: bold;
    }
    @media (max-width: 1200px) {
      .input-field input {
        font-size: 20px !important;
      }
      .sequence-text {
        font-size: 18px;
      }
      .sequence-text span.colored-amino-acid {
        font-size: 18px;
      }
    }
    @media (max-width: 768px) {
      .input-field input {
        font-size: 18px !important;
        height: 55px !important;
      }
      .sequence-text {
        font-size: 16px;
        padding: 12px;
      }
      .sequence-text span.colored-amino-acid {
        font-size: 16px;
        margin-right: 8px;
        min-width: 45px;
      }
      .btn-action {
        width: 100%;
        margin: 0 0 10px 0;
      }
    }
  ")),
  div(class = "intro-text",
      "Welcome, students! 🎓 Observe the generated DNA and RNA sequences. Your task is to enter the correct amino acid sequence by translating the RNA. Let's explore the Central Dogma! 🧬"
  ),
  # Side-by-side layout for input fields (left) and codon chart (right)
  fluidRow(
    # Left column: Input fields and buttons
    column(6, style = "padding-right: 30px; min-width: 320px; max-width: 600px;",
      div(class = "input-section",
        div(class = "input-field",
            tags$label("DNA Template Sequence:"),
            textInput("dna_display", NULL, value = "TACACCGAAGGCTAA")
        ),
        div(class = "input-field",
            tags$label("RNA Sequence:"),
            textInput("rna_display", NULL, value = "AUGUGGCUUCCGAUU")
        ),
        div(class = "input-field",
            tags$label("Enter Amino Acid Sequence (3-letter codes, no spaces):"),
            textInput("amino_acid_input", NULL, value = "", placeholder = "e.g., METTRPLEUPROASPSTOP")
        ),
        div(class = "button-section",
            actionButton("process_sequences", "Check Translation", class = "btn-action btn-process"),
            actionButton("refresh", "Refresh All", class = "btn-action btn-refresh")
        )
      )
    ),
    # Right column: Codon chart and credits
    column(6, style = "padding-left: 30px; min-width: 320px; max-width: 500px;",
      # --- Codon Chart Section ---
      div(class = "chart-section",
        h4("Codon Chart", style = "text-align: center; color: #2c3e50; margin-bottom: 15px;"),
        img(src = "codon_chart.png", style="width: 100%; max-width: 400px; display: block; margin: auto; border-radius: 10px; box-shadow: 0 4px 10px rgba(0,0,0,0.1);"),
        p("Codon chart designed by Dr. Sayumi York & Dr. John Finnerty.", style = "text-align: center; color: #666; font-size: 14px; margin-top: 15px;"),
        p("App coded by Sakib Hussen.", style = "text-align: center; color: #666; font-size: 14px; margin-top: 5px;")
      )
    )
  ),
  # Results section below
  div(class = "results-section",
      h4("Results", style = "color: #2c3e50; text-align: center;"),
      uiOutput("sequence_results"),
      div(class = "error-message", htmlOutput("error_message")),
      div(class = "success-message", htmlOutput("success_message"))
  )
)

# Server Definition
server <- function(input, output, session) {
  # Initial DNA and RNA sequences
  initial_dna <- reactiveVal("TACACCGAAGGCTAA")
  initial_rna <- reactiveVal("AUGUGGCUUCCGAUU") 
  
  # Update display inputs on initial load
  observe({
    updateTextInput(session, "dna_display", value = format_sequence(initial_dna(), chunk_size = 3))
    updateTextInput(session, "rna_display", value = format_sequence(initial_rna(), chunk_size = 3))
  })
  
  # Disable DNA and RNA display inputs to make them readonly
  observe({
    shinyjs::disable("dna_display")
    shinyjs::disable("rna_display")
  })
  
  # Validate inputs and update border colors and show error messages for invalid characters
  observe({
    aa_val <- toupper(gsub(" ", "", input$amino_acid_input))
    aa_invalid_chars <- grepl("[^A-Z]", aa_val) # Check for non-alphabetic for AA
    
    error_msg <- NULL
    if (aa_invalid_chars) {
      error_msg <- "Error: Amino acid sequence contains invalid characters. Only letters are allowed."
    }
    output$error_message <- renderText({
      if (!is.null(error_msg)) {
        error_msg
      } else {
        ""
      }
    })
    
    # Basic input validation for visual feedback
    shinyjs::toggleCssClass("amino_acid_input", "invalid", aa_invalid_chars)
    shinyjs::toggleCssClass("amino_acid_input", "valid", !aa_invalid_chars && nchar(aa_val) > 0)
  })
  
  # Handle translation check
  observeEvent(input$process_sequences, {
    output$error_message <- renderText("")
    output$success_message <- renderText("")
    output$sequence_results <- renderUI(NULL)
    
    current_dna <- initial_dna() 
    current_rna <- initial_rna() 
    
    # First, validate the format of the internally generated RNA (should generally be fine, but good for robustness)
    rna_format_check <- validate_rna_format(current_rna)
    if (grepl("^Error:", rna_format_check)) {
      output$error_message <- renderText(paste("Internal RNA generation error:", rna_format_check))
      return()
    }
    
    # Now, validate the user's amino acid input against the derived RNA
    user_amino_acid_input_processed <- toupper(gsub(" ", "", input$amino_acid_input))
    aa_validation_result <- validate_amino_acids_input(current_rna, user_amino_acid_input_processed)
    
    if (grepl("^Error:", aa_validation_result)) {
      output$error_message <- renderText(aa_validation_result)
      output$sequence_results <- renderUI({
        # Ensure that `translate_rna` is called consistently for correct AA sequence display
        correct_amino_acids_array <- translate_rna(current_rna)
        correct_amino_acids_string <- toupper(paste(correct_amino_acids_array, collapse = ""))
        
        tagList(
          div(class = "sequence-label", "DNA 🧬"),
          div(class = "sequence-text", format_sequence(current_dna, 3)),
          div(class = "sequence-label", "RNA 📜"),
          div(class = "sequence-text", format_sequence(current_rna)),
          div(class = "sequence-label", "Your Amino Acids 🧪"),
          div(class = "sequence-text", format_amino_acids_with_color(user_amino_acid_input_processed)), # Use colored formatting
          div(class = "sequence-label", "Correct Amino Acid Sequence ✅"),
          div(class = "sequence-text", format_amino_acids_with_color(correct_amino_acids_string)) # Use colored formatting
        )
      })
    } else {
      # All sequences are valid and match
      output$error_message <- renderText("")
      output$success_message <- renderText("Congratulations! Your amino acid sequence is correct! 🎉")
      output$sequence_results <- renderUI({
        tagList(
          div(class = "sequence-label", "DNA 🧬"),
          div(class = "sequence-text", format_sequence(current_dna, 3)),
          div(class = "sequence-label", "RNA 📜"),
          div(class = "sequence-text", format_sequence(current_rna)),
          div(class = "sequence-label", "Your Amino Acids 🧪"),
          div(class = "sequence-text", format_amino_acids_with_color(aa_validation_result)) # Use colored formatting
        )
      })
    }
  })
  
  # Handle refresh
  observeEvent(input$refresh, {
    # Always generate a 15-nucleotide DNA sequence for consistency
    new_dna <- generate_dna(length = 15)
    
    # Transcribe the new DNA to get the corresponding RNA
    new_rna <- transcribe_dna(new_dna)
    
    # Double-check that no stop codons are present in the RNA
    if (grepl("UAA|UAG|UGA", new_rna)) {
      # If stop codons found, regenerate the sequence
      for (attempt in 1:10) {
        new_dna <- generate_dna(length = 15)
        new_rna <- transcribe_dna(new_dna)
        if (!grepl("UAA|UAG|UGA", new_rna)) {
          break
        }
      }
      # If still has stop codons, use a completely safe sequence
      if (grepl("UAA|UAG|UGA", new_rna)) {
        safe_bases <- c("A", "C", "G")
        new_dna <- paste(sample(safe_bases, 15, replace = TRUE), collapse = "")
        new_rna <- transcribe_dna(new_dna)
      }
    }
    
    # Update reactive values
    initial_dna(new_dna)
    initial_rna(new_rna)
    
    # Update display inputs
    updateTextInput(session, "dna_display", value = format_sequence(new_dna, chunk_size = 3)) # Format new DNA
    updateTextInput(session, "rna_display", value = format_sequence(new_rna, chunk_size = 3)) # Format new RNA
    updateTextInput(session, "amino_acid_input", value = "")
    output$error_message <- renderText("")
    output$success_message <- renderText("")
    output$sequence_results <- renderUI(NULL)
  })
}

# Run the app
shinyApp(ui = ui, server = server)