"""
A simple BLOSUM matrix reader
"""

class BLOSUM():
    def __init__(self, n=None, path=None):
        """
        Object to easily access a blosum matrix.

        Input:
        Either n ϵ {45,50,62,80,90} or path

        n: Int, which BLOSUM Matrix to use. 
            Choice between: 45,50,62,80 and 90
            Data gathered from https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/
        
        path: String, path to a Blosum matrix.
            File in a format like:
            https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62
        
        """

        self.n = n
        
        if n in [45,50,62,80,90]:
            # Using default matrix
            self.loadMatrix(f"blosum/{n}.blosum")
        elif path is not None:
            # load custom matrix
            self.loadMatrix(path)
        

    def loadMatrix(self, path):
        """
        Reads a Blosum matrix from file.
        File in a format like:
            https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62
        
        Input:
            path: String, path to file
        
        Returns:
            blosumDict: Dictionary, Blosum-Dict
        """
        with open(path, "r") as bFile:
            content = bFile.readlines()

            # Skip header line
            content = content[1:]
            
            # Extract labels
            labels = content[0]
            labels = labels.split()
            
            # Skip label line
            content = content[1:]

            blosumDict = {}
            # For each line
            for line in content:
                line = line.split()
                # Add Line/Label combination to dict
                for index, lab in enumerate(labels, start=1):
                    blosumDict[f"{line[0]}{lab}"] = float(line[index])
        
        self.matrix = blosumDict


    def get(self, aminoA, aminoB):
        """
        Get BLOSUM score
        Return "-inf" if key not found

        Input:
            aminoA: Char
            aminoB: Char
        Ouput:
            score: Float
        """
        score = 0
        try:
            score = self.matrix[f"{aminoA}{aminoB}"]
        except KeyError:
            score = float('-inf')

        return score

