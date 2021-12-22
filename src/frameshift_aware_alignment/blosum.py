"""
A simple BLOSUM matrix reader
"""
from importlib.resources import read_text


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

        # Using default matrix
        if n in [45, 50, 62, 80, 90]:
            path = f"{n}.blosum"
            self.__loadMatrix(path, res=True)

        elif path is not None:
            # load custom matrix
            self.__loadMatrix(path, res=False)
        else:
            Exception("Cant initate empty BLOSUM Object")

    def __loadMatrix(self, path, res=False) -> None:
        """
        Reads a Blosum matrix from file.
        File in a format like:
            https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62

        Input:
            path: str, path to a file.

        Returns:
            blosumDict: Dictionary, Blosum-Dict
        """
        if res:
            content = read_text("frameshift_aware_alignment.blosum_data", path).splitlines()
        else:
            with open(path, "r") as bFile:
                content = bFile.readlines()

        # Skip header line
        content = content[1:]

        # Extract labels
        labels = content[0]
        labelslist = labels.split()

        # Skip label line
        content = content[1:]

        blosumDict = {}
        # For each line
        for line in content:
            linelist = line.split()
            if len(linelist) < 0:
                break
            # Add Line/Label combination to dict
            for index, lab in enumerate(labelslist, start=1):
                blosumDict[f"{linelist[0]}{lab}"] = float(linelist[index])

        self.matrix = blosumDict

    def get(self, aminoA: str, aminoB: str) -> float:
        """
        Get BLOSUM score
        Return "-inf" if key not found

        Input:
            aminoA: Charself.matrix[
            aminoB: Char
        Ouput:
            score: Float
        """
        score = 0.0
        try:
            score = self.matrix[f"{aminoA}{aminoB}"]
        except KeyError:
            score = float('-inf')

        return score
