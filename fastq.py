"""Reading FASTQ files.

Copyright 2022 Maciej Michalczyk (mccv99@gmail.com)

"""

import typing


class Record:
    """Object representing a FASTA (aka Pearson) record.

    Attributes:
        id  (str)         : Sequence identifier
        seq (str)         : Sequence
        description (str) : Description line (defline)
    """

    # saving a whopping 104 bytes from this
    __slots__ = ('id', 'seq', 'desc', 'phred_quality')

    def __init__(self, id: str, seq: str, phread_quality: str, desc: typing.Optional[str] = None):
        """Creates a Record.

        Example:
            >>> record = Record(id='NP_055309.2',
            ...                      seq='MRELEAKAT',
            ...                      desc='TNRC6A')
            >>> print(record)
            >NP_055309.2 TNRC6A
            MRELEAKAT
        """
        self.id = id
        self.seq = seq
        self.desc = desc
        self.phred_quality = phread_quality

    @property
    def description(self) -> str:
        """Returns a description line (defline) of FASTA record.

        Example:
            >>> record = Record(id='NP_055309.2',
            ...                      seq='MRELEAKAT',
            ...                      desc='TNRC6A')
            >>> print(record.description)
            >NP_055309.2 TNRC6A
        """
        lst = [f'>{self.id}']
        if self.desc:
            lst.append(f'{self.desc}')
        return " ".join(lst)

    def __str__(self):
        """Returns the record as a string in the FASTA format.

        Example:
            >>> record = Record(id='NP_055309.2',
            ...                      seq='MRELEAKAT',
            ...                      desc='TNRC6A')
            >>> print(record)
            >NP_055309.2 TNRC6A
            MRELEAKAT
        """
        return self.format(wrap=70).rstrip()

    def format(self, wrap: int = 70):
        """Returns a formatted FASTA record.

        Args:
            wrap:
                Optional line length to wrap sequence lines (default: 70
                characters). Use zero (or None) for no wrapping, giving
                a single long line for the sequence.

        Example:
            >>> record = Record(id='NP_055309.2',
            ...                      seq='MRELEAKAT',
            ...                      desc='TNRC6A')
            >>> print(record.format())
            >NP_055309.2 TNRC6A
            MRELEAKAT
            >>> print(record.format(wrap=3))
            >NP_055309.2 TNRC6A
            MRE
            LEA
            KAT
        """
        lst = [self.description, '\n']
        if wrap:
            for i in range(0, len(self.seq), wrap):
                lst.append(f'{self.seq[i:i + wrap]}\n')
        else:
            lst.append(self.seq)
            lst.append('\n')
            lst.append(self.phred_quality)
        return "".join(lst)
