#Note:BME 160 Final Project
# RT primer Analysis Tool

# This program analyses reverse transcriptase (RT) primers
# It checks how well a primer binds to its intended target gene and how much it might bind to ribosomal RNA (rRNA), which would cause contamination
# The program also suggests alternative primers that maintain strong target binding but reduce rRNA binding

"""
It analyzes reverse transcription (RT) to evaluate potential off-target binding to ribosomal RNA (rRNA).
This tool scans the target gene sequences and rRNA sequences, calculating binding similarity, and suggests alternative
primers that reduce rRNA contamination risk.
"""

"""
def read_fasta(filename):
  - reads fasta file & returns the DNA string as a string

  ex. of fasta file:
  >sequence_name
  ACTGACTGACTG...

  This function will ignore the header line (starting with >), joins all sequence lines together
"""

def read_fasta(filename):
    """ Read fasta file & return the DNA string as a string """
    sequence = ""

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if line == "":
                continue
            if line[0] == ">":
                continue

            # add sequence letters to string
            sequence += line.upper()

    return sequence

#Note:
"""
def reverse_complement(sequence):
  Will return the reverse complement of a DNA sequence.
Primers bind to the opposite DNA strand, so we must search using the reverse complement.

ex. ATCG --> CGAT

"""

def reverse_complement(sequence):
  comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
  new_seq = ""

  #Converts each base to its compliment
  for base in sequence:
    new_seq = new_seq + comp[base]

  #reverse the sequence
  return new_seq[::-1]

#Note:
"""
def count_matches(seq1, seq2):
  - Will count how many bases match between two sequences of the same length
"""

def count_matches(seq1, seq2):
  #counts the number of bases that match between the two sequences of = length
  matches = 0

  for i in range(len(seq1)):
    if seq1[i] == seq2[i]:
      matches += 1

  return matches

#Note:
"""
def percent_match(seq1, seq2):
  - Returns the percent match between two strings
"""

def percent_match(seq1, seq2):
  """calcultes percent match between two sequences"""
  matches = count_matches(seq1, seq2)
  return (matches / len(seq1)) * 100

#Note:
"""
def best_binding_percent(sequence, primer):
  - Return the best percent match of primer anywhere in sequence
Slide the primer across the sequence to find the highest percent match. This simulates checking every possible binding position for the primer.
"""

def best_binding_percent(sequence, primer):

#Slides primer across sequence and finds highest percent match, checks possible binding position for primer
  best_percent = 0

  for i in range(len(sequence) - len(primer) + 1):
    window = sequence[i:i + len(primer)]
    current_percent = percent_match(window, primer)

    if current_percent > best_percent:
      best_percent = current_percent

  return best_percent

#Note:
"""
def generate_mutants(primer):
  - Generates new primers by changing one base at a time.
  Ex.
  Original: ATCG
  mutants: TTCG, CTCG, GTCG, ect.
these are used to test alternative promer design
"""

def generate_mutants(primer):
  #Generates one-base mutant primers
  bases = ["A", "T", "C", "G"]
  mutants = []

  for i in range(len(primer)):
    for base in bases:
      if base != primer[i]:
        new_primer = primer[:i] + base + primer[i +1:]
        mutants.append(new_primer)

  return mutants

#Note:
"""
def suggest_better_primers(original_primer, target_seq, rrna_seq):
  - Will try one base mutant primers and find ones that still bind the target well but reduce rRNA binding
"""

def suggest_better_primers(original_primer, target_seq, rrna_seq):
  #Trys all one base mutant primers with lower rRNA binding, keeps primers that still bind to the target well
  suggestions = []
  mutant_list = generate_mutants(original_primer)

  for primer in mutant_list:
    primer_rc = reverse_complement(primer)

    #target_percent --> checks how well primer binds to target gene
    #rrna_percent --> checks how much it binds to rRNA = contamination risk
    target_percent = best_binding_percent(target_seq, primer_rc)
    rrna_percent = best_binding_percent(rrna_seq, primer_rc)

    #keeps primers that bind to the target well
    if target_percent >= 70:
      suggestions.append([primer, target_percent, rrna_percent])

  #sorts by lowest rRNA first then highest target
  suggestions.sort(key=lambda x: (x[2], -x[1]))

  return suggestions

#Note:
"""
def analyze_primer(primer_name, primer_seq, target_name, target_seq, rrna_seq):
  - Will analyze a single primer and print results

The program will report:
- target binding strength
- rRNA contamination risk
- possible alternative primers
"""

def analyze_primer(primer_name, primer_seq, target_name, target_seq, rrna_seq):
  print()
  print("=" * 60)
  print("Analyzing primer:", primer_name)
  print("Primer sequence:", primer_seq)
  print("=" * 60)

  #convert primer to reverse complement
  primer_rc = reverse_complement(primer_seq)
  print("Reverse complement used for search:", primer_rc)

  #Calculate best matches
  target_percent = best_binding_percent(target_seq, primer_rc)
  rrna_percent = best_binding_percent(rrna_seq, primer_rc)

  print("Best on-target sequence match to", target_name, "=", round(target_percent, 2), "%")
  print("Best sequence match to rRNA (off-target risk) =", round(rrna_percent, 2), "%")

  #determines if primer passes contamination threshold
  if rrna_percent <= 30:
    print("This primer passes the off-target rule.")
  else:
    print("This primer fails the off-target rule.")

  print()
  print("Suggested alternative primers:")

  suggestions = suggest_better_primers(primer_seq, target_seq, rrna_seq)

  if len(suggestions) == 0:
    print("No alternative primers found.")
  else:
    for suggestion in suggestions [:10]:
      print("Primer:", suggestion[0], "Target match:", round(suggestion[1], 2), "%", "rRNA Match:", round(suggestion[2], 2), "%")

    print()
    best_suggestion = suggestions[0]
    print("Best alternative found:")
    print("Primer:", best_suggestion[0])
    print("Target match:", round(best_suggestion[1], 2), "%")
    print("rRNA match = ", round(best_suggestion[2], 2), "%")
    print("GC content: " + str(GC_Content(best_suggestion[0])) + "% GC")

    if best_suggestion[2] <= 30:
       print("This suggestion passes the 30% off-target rule.")
    else:
       print("This suggestion improves the primer but fails the 30% off-target rule.")

def GC_Content(in_str):
  gc = 0
  for i in in_str:
    if i == "G" or i == "C":
      gc +=1

  return gc/len(in_str)*100

#Note:
"""
def main():
  - Main fucntion that loads sequence and analyzes primers
"""

def main():
  # read fasta sequence files
  xbp1_seq = read_fasta("xbp1.fasta")
  ets4_seq = read_fasta("ets4.fasta")
  rrna_seq = read_fasta("rrna.fasta")

  # prints sequence lengths to confirm files loaded correctly
  print("Length of xbp1 sequence =", len(xbp1_seq))
  print("Length of ets4 sequence =", len(ets4_seq))
  print("Length of rrna sequence =", len(rrna_seq))
  print()

  # primers form the experimental dataset
  xbp1_primer = "CCTCAATGTAGCCAGCAACG"
  ets4_primer = "AGGCTAGTAAACAGGGGAAGGA"

  # prints reverse complements for verification
  print("xbp-1 reverse complement =", reverse_complement(xbp1_primer))
  print("ets-4 reverse complement =", reverse_complement(ets4_primer))
  print()

  # checks that the target sequnce contains primer binding sites
  print("Does xbp1 contain xbp1 primer reverse compliment?", reverse_complement(xbp1_primer) in xbp1_seq)
  print("Does ets4 contain ets4 primer reverse compliment?", reverse_complement(ets4_primer) in ets4_seq)
  print()

  # analyze each primer
  analyze_primer("xbp-1 primer", xbp1_primer, "xbp-1", xbp1_seq, rrna_seq)
  analyze_primer("ets-4 primer", ets4_primer, "ets-4", ets4_seq, rrna_seq)


# runs the program
main()