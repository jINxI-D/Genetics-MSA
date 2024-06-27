import sys
import csv
import os
import subprocess
import matplotlib.pyplot as plt
import tempfile
from scipy.stats import hypergeom
from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QFileDialog,
    QMessageBox,
    QTableWidget,
    QTableWidgetItem,
    QCheckBox,
)
from PyQt5.QtCore import Qt


class ConservationAnalyzer(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Conservation Analyzer")

        # Input fields
        self.alignment_file_label = QLabel("Alignment File:")
        self.alignment_file_input = QLineEdit()
        self.alignment_file_button = QPushButton("Browse")
        self.alignment_file_button.clicked.connect(self.open_alignment_file)

        self.conservation_threshold_label = QLabel("Conservation Threshold (0-1):")
        self.conservation_threshold_input = QLineEdit("0.8")

        self.mutation_threshold_label = QLabel("Mutation Threshold (0-1):")
        self.mutation_threshold_input = QLineEdit("0.2")

        self.aligner_label = QLabel("Aligner:")
        self.aligner_mafft_checkbox = QCheckBox("MAFFT")
        self.aligner_mafft_checkbox.setChecked(True)
        self.aligner_clustalo_checkbox = QCheckBox("Clustal Omega")

        # Reference sequence input
        self.reference_sequence_label = QLabel("Reference Sequence ID (optional):")
        self.reference_sequence_input = QLineEdit()

        # Output section
        self.output_label = QLabel("Output:")
        self.conserved_sites_table = QTableWidget()
        self.conserved_sites_table.setColumnCount(3)
        self.conserved_sites_table.setHorizontalHeaderLabels(["Position", "Residue", "Conservation Rate"])
        self.conserved_sites_table.verticalHeader().setVisible(False)

        self.mutated_sites_table = QTableWidget()
        self.mutated_sites_table.setColumnCount(3)
        self.mutated_sites_table.setHorizontalHeaderLabels(["Position", "Residue", "Conservation Rate"])
        self.mutated_sites_table.verticalHeader().setVisible(False)

        # Analysis button
        self.analyze_button = QPushButton("Analyze")
        self.analyze_button.clicked.connect(self.analyze)

        # Layout
        main_layout = QVBoxLayout()

        # Input section
        input_layout = QVBoxLayout()
        input_layout.addWidget(self.alignment_file_label)
        input_layout.addWidget(self.alignment_file_input)
        input_layout.addWidget(self.alignment_file_button)

        input_layout.addWidget(self.conservation_threshold_label)
        input_layout.addWidget(self.conservation_threshold_input)

        input_layout.addWidget(self.mutation_threshold_label)
        input_layout.addWidget(self.mutation_threshold_input)

        aligner_layout = QHBoxLayout()
        aligner_layout.addWidget(self.aligner_label)
        aligner_layout.addWidget(self.aligner_mafft_checkbox)
        aligner_layout.addWidget(self.aligner_clustalo_checkbox)
        input_layout.addLayout(aligner_layout)

        input_layout.addWidget(self.reference_sequence_label)
        input_layout.addWidget(self.reference_sequence_input)

        main_layout.addLayout(input_layout)

        # Output section
        output_layout = QVBoxLayout()
        output_layout.addWidget(self.output_label)
        output_layout.addWidget(self.conserved_sites_table)
        output_layout.addWidget(self.mutated_sites_table)
        main_layout.addLayout(output_layout)

        # Analysis button
        main_layout.addWidget(self.analyze_button)

        self.setLayout(main_layout)

    def open_alignment_file(self):
        """Opens a file dialog to select the alignment file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Alignment File", "", "FASTA Files (*.fasta)"
        )
        if file_path:
            self.alignment_file_input.setText(file_path)

    def analyze(self):
    # Aligner choice validation
        if self.aligner_mafft_checkbox.isChecked() == self.aligner_clustalo_checkbox.isChecked():
            QMessageBox.critical(self, "Error", "Please select exactly one aligner.")
            return

        alignment_file = self.alignment_file_input.text()
        conservation_threshold = float(self.conservation_threshold_input.text())
        mutation_threshold = float(self.mutation_threshold_input.text())
        aligner = "mafft" if self.aligner_mafft_checkbox.isChecked() else "clustalo"
        reference_sequence_id = self.reference_sequence_input.text()
        try:

            if not (
                0 <= conservation_threshold <= 1
                and 0 <= mutation_threshold <= 1
            ):
                raise ValueError("Thresholds must be between 0 and 1")

            sequences = self.read_fasta(alignment_file)

            if not self.is_aligned(sequences):
                print("Input sequences are not aligned. Performing multiple sequence alignment...")
                aligned_sequences = self.perform_msa(alignment_file, aligner)
                with open("aligned.fasta", "w") as aligned_file:
                    aligned_file.write(aligned_sequences)
                sequences = self.read_fasta("aligned.fasta")

            overall_rate, position_rates = self.calculate_conservation_rate(
                sequences, ref_seq_index=self.find_reference_index(sequences, reference_sequence_id)
            )

            conserved_positions = self.filter_by_rate(
                position_rates, conservation_threshold, is_conserved=True
            )
            mutated_positions = self.filter_by_rate(
                position_rates, mutation_threshold, is_conserved=False
            )

            # Statistical significance test (hypergeometric test)
            num_sites = len(position_rates)
            num_conserved = len(conserved_positions)
            num_mutated = len(mutated_positions)
            p_value_conserved = self.calculate_p_value(
                num_sites, num_conserved, num_mutated
            )
            p_value_mutated = self.calculate_p_value(
                num_sites, num_mutated, num_conserved
            )

            # Update GUI tables
            self.update_table(self.conserved_sites_table, conserved_positions)
            self.update_table(self.mutated_sites_table, mutated_positions)

            # Display results
            QMessageBox.information(
                self,
                "Analysis Results",
                f"Overall conservation rate: {overall_rate:.4f}\n"
                f"Number of conserved positions: {len(conserved_positions)} (p-value: {p_value_conserved:.4f})\n"
                f"Number of highly mutated positions: {len(mutated_positions)} (p-value: {p_value_mutated:.4f})\n"
            )

            # Plot conservation rates
            self.plot_conservation_rates(
                position_rates, conservation_threshold, mutation_threshold
            )

        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def read_fasta(self, file_path):
        """Reads FASTA sequences from a file."""
        sequences = []
        with open(file_path, "r") as file:
            current_id = None
            current_seq = ""
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if current_id is not None:
                        sequences.append((current_id, current_seq))
                    current_id = line[1:]
                    current_seq = ""
                else:
                    current_seq += line
            if current_id is not None:
                sequences.append((current_id, current_seq))
        return sequences

    def is_aligned(self, sequences):
        """Checks if sequences are aligned."""
        if not sequences:
            return True
        first_length = len(sequences[0][1])
        return all(len(seq[1]) == first_length for seq in sequences)

    def perform_msa(self, input_file, aligner="mafft"):
        aligned_sequences = ""  # Initialize the variable to ensure it has a scope beyond the try-except block
        try:
            with tempfile.NamedTemporaryFile(mode='r+', delete=False) as temp_aligned_file:
                if aligner == "mafft":
                    subprocess.run(["mafft", "--quiet", input_file], stdout=temp_aligned_file, text=True, check=True)
                elif aligner == "clustalo":
                    subprocess.run(["clustalo", "-i", input_file, "--outfmt=fasta", "--force"], stdout=temp_aligned_file, text=True, check=True)
                else:
                    raise ValueError("Invalid aligner specified. Choose either 'mafft' or 'clustalo'.")
                temp_aligned_file_path = temp_aligned_file.name
            # Read from the temporary file after the subprocess completes and the file is closed
            with open(temp_aligned_file_path, 'r') as file:
                aligned_sequences = file.read()
        except FileNotFoundError as e:
            print(f"Error: {aligner.upper()} not found. Please install {aligner.upper()} or choose another aligner.")
            raise e  # or handle it as needed
        except subprocess.CalledProcessError as e:
            print(f"Error: Alignment with {aligner.upper()} failed.")
            raise e  # or handle it as needed
        finally:
            # Clean up the temporary file if it exists
            if os.path.exists(temp_aligned_file_path):
                os.remove(temp_aligned_file_path)
        return aligned_sequences

    def calculate_conservation_rate(self, sequences, ref_seq_index=None):
        """Calculates conservation rates."""
        if not sequences:
            return 0, []

        alignment_length = len(sequences[0][1])
        num_sequences = len(sequences)
        conservation_rates = []

        for i in range(alignment_length):
            column = [seq[i] for _, seq in sequences]

            if ref_seq_index is not None:
                reference_residue = sequences[ref_seq_index][1][i]
            else:
                reference_residue = max(set(column), key=column.count)

            conservation_rate = round(
                column.count(reference_residue) / num_sequences, 4
            )
            conservation_rates.append(conservation_rate)

        overall_conservation_rate = round(
            sum(conservation_rates) / alignment_length, 4
        )
        return overall_conservation_rate, conservation_rates

    def filter_by_rate(self, position_rates, threshold, is_conserved=True):
        """Filters positions based on conservation rate."""
        filtered_positions = []
        for i, rate in enumerate(position_rates):
            if (is_conserved and rate >= threshold) or (
                not is_conserved and rate <= threshold
            ):
                filtered_positions.append((i + 1, rate))
        return filtered_positions

    def find_reference_index(self, sequences, reference_sequence_id):
        """Finds the index of the reference sequence."""
        for i, (seq_id, _) in enumerate(sequences):
            if seq_id == reference_sequence_id:
                return i
        return None

    def calculate_p_value(self, num_sites, num_successes, num_failures):
        """Calculates the p-value using the hypergeometric test."""
        return hypergeom.sf(num_successes - 1, num_sites, num_successes, num_failures)

    def get_residue_at_position(self, position):
        sequences = self.read_fasta(self.alignment_file_input.text())
        alignment_length = len(sequences[0][1]) if sequences else 0

        if position <= 0 or position > alignment_length:
            return "N/A"  # Position out of range for the alignment

        column = [seq[position - 1] for _, seq in sequences if len(seq) >= position]

        # Handling when no reference sequence ID is specified or it doesn't match
        if not self.reference_sequence_input.text().strip() or all(seq_id != self.reference_sequence_input.text().strip() for seq_id, _ in sequences):
            # Find the most frequent residue, excluding '-'
            most_frequent_residue = max(set(column) - {"-"}, key=column.count, default="N/A")
            return most_frequent_residue
        else:
            # If a specific reference sequence ID is provided and found
            for seq_id, seq in sequences:
                if seq_id == self.reference_sequence_input.text().strip() and len(seq) >= position:
                    return seq[position - 1]  # Return the residue at the specified position

        return "N/A"  # Fallback if the residue cannot be determined

    def update_table(self, table, data):
        """Updates a QTableWidget with data."""
        table.setRowCount(len(data))
        for row, (position, rate) in enumerate(data):
            table.setItem(row, 0, QTableWidgetItem(str(position)))
            residue = self.get_residue_at_position(position)
            if residue != "N/A":
                table.setItem(row, 1, QTableWidgetItem(str(residue)))
                table.setItem(row, 2, QTableWidgetItem(f"{rate:.4f}"))

    def plot_conservation_rates(
        self, position_rates, conservation_threshold, mutation_threshold
    ):
        """Plots conservation rates."""
        plt.bar(range(1, len(position_rates) + 1), position_rates)
        plt.xlabel("Position")
        plt.ylabel("Conservation Rate")
        plt.title("Conservation Rates")
        plt.axhline(
            y=conservation_threshold, color="g", linestyle="--", label="Conservation Threshold"
        )
        plt.axhline(
            y=mutation_threshold, color="r", linestyle="--", label="Mutation Threshold"
        )
        plt.legend()
        plt.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    analyzer = ConservationAnalyzer()
    analyzer.show()
    sys.exit(app.exec_())
