
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
  "ALA" = list(category = "Aliphatic", color = "#FF0000"), # Alanine
  "GLY" = list(category = "Aliphatic", color = "#FF0000"), # Glycine
  "ILE" = list(category = "Aliphatic", color = "#FF0000"), # Isoleucine
  "LEU" = list(category = "Aliphatic", color = "#FF0000"), # Leucine
  "PRO" = list(category = "Aliphatic", color = "#FF0000"), # Proline
  "VAL" = list(category = "Aliphatic", color = "#FF0000"), # Valine
  
  # Aromatic (Green) - Adjusted for better visibility
  "PHE" = list(category = "Aromatic", color = "#228B22"), # Phenylalanine (Forest Green)
  "TRP" = list(category = "Aromatic", color = "#228B22"), # Tryptophan
  "TYR" = list(category = "Aromatic", color = "#228B22"), # Tyrosine
  
  # Acidic (Blue)
  "ASP" = list(category = "Acidic", color = "#0000FF"), # Aspartic Acid
  "GLU" = list(category = "Acidic", color = "#0000FF"), # Glutamic Acid
  
  # Basic (Light Blue/Cyan) - Adjusted for better visibility
  "ARG" = list(category = "Basic", color = "#4682B4"), # Arginine (Steel Blue)
  "HIS" = list(category = "Basic", color = "#4682B4"), # Histidine
  "LYS" = list(category = "Basic", color = "#4682B4"), # Lysine
  
  # Hydroxylic (Pink)
  "SER" = list(category = "Hydroxylic", color = "#FF69B4"), # Serine (Hot Pink)
  "THR" = list(category = "Hydroxylic", color = "#FF69B4"), # Threonine
  
  # Sulfur-containing (Orange)
  "CYS" = list(category = "Sulfur-containing", color = "#FFA500"), # Cysteine
  "MET" = list(category = "Sulfur-containing", color = "#FFA500"), # Methionine
  
  # Amidic (Purple)
  "ASN" = list(category = "Amidic", color = "#800080"), # Asparagine
  "GLN" = list(category = "Amidic", color = "#800080"), # Glutamine
  
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

# Color mapping for amino acids (3-letter codes) based on properties
amino_acid_properties <- list(
  "ALA" = list(category = "Aliphatic", color = "#FF0000"), # Alanine - Red
  "GLY" = list(category = "Aliphatic", color = "#FF0000"), # Glycine - Red
  "ILE" = list(category = "Aliphatic", color = "#FF0000"), # Isoleucine - Red
  "LEU" = list(category = "Aliphatic", color = "#FF0000"), # Leucine - Red
  "PRO" = list(category = "Aliphatic", color = "#FF0000"), # Proline - Red
  "VAL" = list(category = "Aliphatic", color = "#FF0000"), # Valine - Red
  "PHE" = list(category = "Aromatic", color = "#228B22"),  # Phenylalanine - Green
  "TRP" = list(category = "Aromatic", color = "#228B22"),  # Tryptophan - Green
  "TYR" = list(category = "Aromatic", color = "#228B22"),  # Tyrosine - Green
  "ASP" = list(category = "Acidic", color = "#0000FF"),    # Aspartic Acid - Blue
  "GLU" = list(category = "Acidic", color = "#0000FF"),    # Glutamic Acid - Blue
  "ARG" = list(category = "Basic", color = "#4682B4"),     # Arginine - Steel Blue
  "HIS" = list(category = "Basic", color = "#4682B4"),     # Histidine - Steel Blue
  "LYS" = list(category = "Basic", color = "#4682B4"),     # Lysine - Steel Blue
  "SER" = list(category = "Hydroxylic", color = "#FF69B4"),# Serine - Hot Pink
  "THR" = list(category = "Hydroxylic", color = "#FF69B4"),# Threonine - Hot Pink
  "CYS" = list(category = "Sulfur-containing", color = "#FFA500"), # Cysteine - Orange
  "MET" = list(category = "Sulfur-containing", color = "#FFA500"), # Methionine - Orange
  "ASN" = list(category = "Amidic", color = "#800080"),    # Asparagine - Purple
  "GLN" = list(category = "Amidic", color = "#800080"),    # Glutamine - Purple
  "STOP" = list(category = "Stop Codon", color = "#696969")# Stop Codon - Dim Gray
)

# Function to format amino acid sequence with colors based on properties
format_amino_acids_with_color <- function(amino_acid_string) {
  if (is.null(amino_acid_string) || nchar(amino_acid_string) == 0) {
    return("")
  }
  
  amino_acid_string <- toupper(gsub(" ", "", amino_acid_string))
  # Split into 3-letter chunks
  amino_acids <- strsplit(amino_acid_string, "(?<=\\G...)", perl = TRUE)[[1]]
  
  # Generate HTML for each amino acid with its corresponding color
  formatted_html_elements <- lapply(amino_acids, function(aa_code) {
    props <- amino_acid_properties[[aa_code]]
    if (is.null(props)) {
      # Fallback for unrecognized amino acids or if properties are not defined
      return(tags$span(style = "color: #333333;", aa_code)) # Default to dark gray
    } else {
      return(tags$span(style = paste0("color: ", props$color, "; font-weight: bold; margin-right: 8px;"), aa_code))
    }
  })
  
  # Combine the individual span elements with a non-breaking space
  # This ensures spaces are rendered correctly between colored chunks
  do.call(tagList, lapply(formatted_html_elements, function(x) list(x, tags$span("&nbsp;&nbsp;"))))
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
  paste(substring(seq_string, seq(1, nchar(seq_string), chunk_size), seq(chunk_size, nchar(seq_string), chunk_size)), collapse = " ")
}

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
      return(tags$span(class = "colored-amino-acid", style = "color: #333333;", aa_code)) # Add a class
    } else {
      return(tags$span(class = "colored-amino-acid", style = paste0("color: ", props$color, ";"), aa_code)) # Add a class
    }
  })
  
  # Simply combine the individual span elements. CSS will handle spacing.
  # Use HTML() to ensure the list of tags is rendered as HTML, rather than R's default list printing.
  # You might need to flatten the list if it's deeply nested, but `tagList` usually handles this.
  tagList(formatted_html_elements)
}


# Function to generate random DNA sequence
generate_dna <- function(length = 15) {
  bases <- c("A", "T", "C", "G")
  paste(sample(bases, length, replace = TRUE), collapse = "")
}

# UI Definition
ui <- fluidPage(
  useShinyjs(),
  titlePanel("🌟 Central Dogma Explorer"),
  tags$style(HTML("
    body {
      background-color: #f0f4f8;
      font-family: 'Inter', sans-serif; /* Changed to Inter */
      color: #333;
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
    .shiny-input-container {
      margin-bottom: 20px;
    }
    input[type='text'] {
      font-family: 'Courier New', monospace !important;
      font-size: 24px !important;
      /* letter-spacing: 3px !important; */ /* Removed to prevent interference with colored spans */
      padding: 10px !important;
      border-radius: 8px !important; /* Rounded corners */
      border: 2px solid #ccc !important;
      width: 100% !important; /* Made responsive */
      max-width: 600px; /* Max width for larger screens */
      height: 50px !important;
      box-sizing: border-box;
      overflow-x: auto !important;
      white-space: nowrap !important;
    }
    /* Color-code nucleotides */
    input[type='text']::placeholder {
      color: #999;
    }
    .valid {
      border-color: #28a745 !important;
    }
    .invalid {
      border-color: #dc3545 !important;
    }
    .btn-action {
      font-size: 16px;
      padding: 10px 20px;
      border-radius: 5px;
      margin-right: 10px;
      transition: background-color 0.3s ease-in-out, transform 0.2s ease-in-out; /* Smooth transitions */
      box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    .btn-action:active {
        transform: translateY(1px); /* Slight press effect */
    }
    .btn-process {
      background-color: #007bff;
      color: white;
      border: none;
    }
    .btn-process:hover {
      background-color: #0056b3;
    }
    .btn-refresh {
      background-color: #28a745;
      color: white;
      border: none;
    }
    .btn-refresh:hover {
      background-color: #218838;
    }
    .sequence-label {
      font-weight: bold;
      font-size: 18px;
      color: #2c3e50;
      margin-top: 15px;
    }
    .sequence-text {
      font-family: 'Courier New', monospace;
      font-size: 20px;
      etter-spacing: normal; /* Ensure default letter spacing for the container */
      white-space: normal; /* Allow wrapping if sequence is too long for one line */
      word-break: break-all; 
      background-color: #fff;
      padding: 10px;
      border-radius: 5px;
      box-shadow: 0 2px 5px rgba(0,0,0,0.1);
      word-wrap: break-word; 
      overflow-x: auto;
      
      max-width: 100%;
      min-height: 40px; /* Ensure a minimum height even if empty */
    }
      .sequence-text span.colored-amino-acid {
      margin-right: 15px; /* Adjust this value for more or less space between codes */
      font-family: 'Courier New', monospace; /* Ensure the font is applied to the span */
      letter-spacing: normal; /* Reset letter-spacing if it was inherited differently */
      display: inline-block; /* Allows margin-right to work correctly */
      min-width: 45px; /* Ensure consistent width for each 3-letter code */
      text-align: left;
    }
    .error-message {
      color: #dc3545;
      font-size: 16px;
      margin-top: 10px;
      font-weight: bold;
    }
    .success-message {
      color: #28a745;
      font-size: 16px;
      margin-top: 10px;
      font-weight: bold;
    }
    .sidebar-panel {
      background-color: #ffffff;
      padding: 20px;
      border-radius: 10px;
      box-shadow: 0 4px 10px rgba(0,0,0,0.1);
      margin-bottom: 20px; /* Add some spacing */
    }
    .main-panel {
      background-color: #ffffff;
      padding: 20px;
      border-radius: 10px;
      box-shadow: 0 4px 10px rgba(0,0,0,0.1);
    }
    /* Responsive adjustments */
    @media (max-width: 768px) {
      input[type='text'] {
        width: 100% !important;
        max-width: none;
        font-size: 18px !important;
      }
      .sequence-text {
        font-size: 16px;
      }
      .btn-action {
        width: 100%;
        margin-right: 0;
        margin-bottom: 10px;
      }
    }
  ")),
  div(class = "intro-text",
      "Welcome, students! 🎓 Observe the generated DNA and RNA sequences. Your task is to enter the correct amino acid sequence by translating the RNA. Let’s explore the Central Dogma! 🧬"
  ),
  div(class = "container-fluid", # Use container-fluid for better responsiveness
      fluidRow(
        column(12, class = "sidebar-panel",
               textInput("dna_display", "DNA Template Sequence:", value = "TACACCGAAGGCTAA"),
               textInput("rna_display", "RNA Sequence:", value = "AUGUGGCUUCCGAUU"),
               textInput("amino_acid_input", "Enter Amino Acid Sequence (3-letter codes, no spaces):", value = "",
                         placeholder = "e.g., METTRPLEUPROASPSTOP"),
               div(style = "text-align: center;", # Center buttons
                   actionButton("process_sequences", "Check Translation", class = "btn-action btn-process"),
                   actionButton("refresh", "Refresh All", class = "btn-action btn-refresh")
               )
        )
      ),
      fluidRow(
        column(12, class = "main-panel",
               style = "margin-top: 20px;",
               h4("Results", style = "color: #2c3e50; text-align: center;"),
               uiOutput("sequence_results"),
               div(class = "error-message", htmlOutput("error_message")),
               div(class = "success-message", htmlOutput("success_message"))
        )
      )
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
          div(class = "sequence-text", format_sequence(current_dna, 1)),
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
          div(class = "sequence-text", format_sequence(current_dna, 1)),
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
    # Generate a new DNA sequence (length 15-21, multiple of 3)
    new_dna_length <- sample(seq(15, 21, by = 3), 1)
    new_dna <- generate_dna(length = new_dna_length)
    
    # Transcribe the new DNA to get the corresponding RNA
    new_rna <- transcribe_dna(new_dna)
    
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
