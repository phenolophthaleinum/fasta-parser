"""Reading and writing FASTA files.

Copyright 2022 Andrzej Zielezinski (a.zielezinski@gmail.com)
https://github.com/aziele/fasta-parser
"""

import bz2
import gzip
import pathlib
import typing
import zstandard
import lz4.frame as lz4
import mmap


class Record:
    """Object representing a FASTA (aka Pearson) record.

    Attributes:
        id  (str)         : Sequence identifier
        seq (str)         : Sequence
        description (str) : Description line (defline)
    """

    # saving a whopping 104 bytes from this
    __slots__ = ('id', 'seq', 'desc')

    def __init__(self, id: str, seq: str, desc: typing.Optional[str] = None):
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

    def __iter__(self):
        """Iterates over the letters in the sequence.

        Example:
            >>> record = Record(id='NP_055309.2', 
            ...                      seq='MRELEAKAT',
            ...                      desc='TNRC6A')
            >>> for amino_acid in record:
            ...     print(amino_acid)
            M
            R
            E
            L
            E

            This is equivalent to iterating over the sequence directly:
            >>> for amino_acid in record.seq:
            ...     print(amino_acid)
            M
            R
            E
            L
            E
        """
        return iter(self.seq)

    def __contains__(self, char):
        """Implements the 'in' keyword to search the sequence.

        Example:
            >>> record = Record(id='NP_055309.2', 
            ...                      seq='MRELEAKAT',
            ...                      desc='TNRC6A')
            >>> print('M' in record)
            True
        """
        return char in self.seq

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

    def __len__(self):
        """Return the length of the sequence.

        Example:
            >>> record = Record(id='NP_055309.2',
            ...                      seq='MRELEAKAT',
            ...                      desc='TNRC6A')
            >>> len(record)
            9
        """
        return len(self.seq)

    def format(self, wrap:int = 70):
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
        return "".join(lst)


def parse(filename: typing.Union[str, pathlib.Path]):
    """Iterates over FASTA records in a file.

    Args:
        filename: A name or path of file containing FASTA sequences.

    Returns:
        A generator of Record objects.
    """
    seqid = None
    desc = None
    seq = []
    with get_open_func(filename)(filename, 'rt') as fh:
        # with mmap.mmap(fh.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_fh:
        #     print(mmap_fh[:50])
        for line in fh:
            if line.startswith('>'):
                if seq:
                    yield Record(seqid, "".join(seq), desc)
                    seq = []
                seqid = line.split()[0][1:]
                desc = line[len(seqid)+1:].strip()
            else:
                seq.append(line.strip())
        if seq:
            yield Record(seqid, "".join(seq), desc)


def to_dict(sequences):
    """Turns a generator or list of Record objects into a dictionary.

    This function is not suitable for very large sets of sequences, as all the
    SeqRecord objects are held in memory.

    Args:
        sequences: an iterator that returns Record objects, or simply a 
          list of SeqRecord objects.

    Returns:
        A dict mapping sequence id (key) to Record object (value).

    Raises:
        If there are duplicate keys, an error is raised.

    Example:
        >>> import Fasta
        >>> pdict = Fasta.to_dict(Fasta.parse('test.fa'))
        >>> print(sorted(pdict.keys()))
        ['gi|195354411|', 'tr|Q8SY33|']
        >>> print(pdict['tr|Q8SY33|'].description)
        Gawky, isoform A [Drosophila melanogaster]
        >>> len(pdict)
        2
    """
    return {record.id: record for record in sequences}


def get_compression_type(filename: typing.Union[str, pathlib.Path]) -> str:
    """Guesses the compression (if any) on a file using the first few bytes.

    http://stackoverflow.com/questions/13044562

    Returns:
        Compression type (gz, bz2, plain)
    """
    magic_dict = {(b'\x1f', b'\x8b', b'\x08'): 'gz',
                  (b'\x42', b'\x5a', b'\x68'): 'bz2',
                  (b'\x50', b'\x4b', b'\x03', b'\x04'): 'zip',
                  b'(\xb5/\xfd': 'zst',
                  b'\x04"M\x18': 'lz4'}
    # since recognized compressions are not added dynamically, the max size of magic bytes can be static
    max_len = 4

    fh = open(str(filename), 'rb')
    file_start = fh.read(max_len)
    fh.close()
    compression_type = [magic_dict[elem] for elem in magic_dict if file_start.startswith(elem)]

    return compression_type[0] if compression_type else 'plain'


def get_open_func(filename: typing.Union[str, pathlib.Path]):
    """Returns function to open a file."""
    open_funcs = {
        'gz': gzip.open,
        'bz2': bz2.open,
        'plain': open,
        'zst': zstandard.open,
        'lz4': lz4.open
    }

    return open_funcs[get_compression_type(filename)]
