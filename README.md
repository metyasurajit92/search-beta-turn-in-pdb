This program is designed to identify and analyze single and double beta turns within a given PDB (Protein Data Bank) structure. Beta turns are a key secondary structural feature in proteins, often involved in reversing the direction of the polypeptide chain and contributing to the compact, globular nature of proteins. 

By utilizing the BioStructure.jl library, this program provides a robust and efficient framework for parsing and interpreting the complex three-dimensional structural data contained in PDB files. The BioStructure.jl library is a powerful tool for structural bioinformatics in Julia, offering functionalities for reading PDB files, manipulating atomic coordinates, and identifying structural motifs.

Key Features:
1. Automated Detection of Beta Turns: 
   - Identifies regions within the protein structure that correspond to single beta turns (one loop involving four residues) and double beta turns (two consecutive turns).
   
2. Structural Analysis: 
   - Calculates geometric parameters such as torsion angles (phi, psi) and interatomic distances to confirm the presence of beta turns.
   - Provides insights into the conformation of the identified turns and their roles in the overall protein structure.

3. User-Friendly Input and Output: 
   - Accepts a standard PDB file as input, allowing seamless integration with existing protein structural data.
   - Outputs detailed reports on the identified beta turns, including residue indices, sequence information, and geometric properties.

4. Scalability and Flexibility: 
   - Capable of analyzing large and complex protein structures, making it suitable for both small-scale academic studies and high-throughput bioinformatics workflows.

### Implementation Details:
- The program leverages the parsing capabilities of BioStructure.jl to efficiently extract atomic and residue-level data from PDB files.
- Using Julia's high-performance computation capabilities, it processes structural data with minimal computational overhead.
- The beta turn detection algorithm is implemented to adhere to established criteria, such as the Ramachandran plot values and specific inter-residue hydrogen bonding patterns.

This tool can be invaluable for structural biologists and bioinformaticians interested in protein folding, stability, and function, as beta turns often play a critical role in protein interactions and dynamics.
