"""
@author: Clara Heimburger, Jules Kreuer
"""

from os import path


class fasta_object():
    def __init__(self, head, body):
        """
        Object to keep a valid fasta object
        """
        self.head = head
        self.body = body


class fasta_handler():
    def read(self, fileName):
        """
        Reads a fasta-style file and returns a list of fasta_objects
        """

        if not path.isfile(fileName):
            raise FileNotFoundError("FastaFile not found!")

        fasta_objects = []
        with open(fileName, 'r') as f:
            head = ""
            body = ""
            newObject = True
            for line in f:
                # go through each line
                # First Header
                if newObject and line.startswith(">"):
                    head = line.strip()
                    body = ""
                    newObject = False
                # N-th Header
                elif line.startswith(">"):
                    fasta_objects.append(fasta_object(head, body))
                    head = line.strip()
                    body = ""
                # Sequence
                else:
                    body += line.strip().upper()
            # append last element
            fasta_objects.append(fasta_object(head, body))
        return fasta_objects

    def write(self, fasta_pairs, fileName, mode="w"):
        """
        Writes a list of fasta_objects or a single one to a file.
        Takes fasta_objects as input.
        """

        if not isinstance(fasta_pairs, list):
            fasta_pairs = [fasta_pairs]

        with open(fileName, mode) as f:
            for fo in fasta_pairs:
                f.write(fo.head + "\n")
                body_len = len(fo.body)
                # Write only 70 chars per line
                for i in range(0, body_len, 70):
                    f.write(fo.body[i:i+70] + "\n")
        return None


    def print(self, fasta):
        """
        Prints a single or a list of fasta_objects.
        Takes fasta_objects as input.
        """

        if not isinstance(fasta, list):
            fasta = [fasta]
        
        for fo in fasta:
            print(fo.head)
            body_len = len(fo.body)
            # Print only 70 chars per line
            for i in range(0, body_len, 70):
                print(fo.body[i:i+70])
        return None
             
            
