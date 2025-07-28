from sys import exit
from typing import Any, Generator, Optional
from functools import cache

all = ["Minigfa", "Minibed"]


class _Segment:
    def __init__(self, S_line: str) -> None:
        easy_line = S_line.strip().split("\t")
        self.segID: str = easy_line[1]
        self.seq: str = easy_line[2]
        SN_parts = easy_line[4].split(":")[2].split("#")
        self._SName = (SN_parts[0], int(SN_parts[1]), SN_parts[2])
        self.SRank: int = int(easy_line[6].split(":")[2])

    @property
    def source_sample(self) -> str:
        return self._SName[0]

    @property
    def linear_reference_chr(self) -> str:
        return self._SName[2]


class _Link:
    def __init__(self, L_line: str) -> None:
        easy_line = L_line.strip().split("\t")
        self.fromSeg = (easy_line[1], easy_line[2])
        self.toSeg = (easy_line[3], easy_line[4])
        self.SRank = int(easy_line[6].split(":")[2])

    @property
    def fromID(self) -> str:
        return self.fromSeg[0]

    @property
    def fromOrient(self) -> str:
        return self.fromSeg[1]

    @property
    def toID(self) -> str:
        return self.toSeg[0]

    @property
    def toOrient(self) -> str:
        return self.toSeg[1]


class Minigfa:
    """
        Composite data storing GFA file information.
        Notice:Only lines S and L can be processed, lines starting with other letters are discarded

        The path of the GFA file that is preferably passed in when constructing the object.If you don't do this, you will just get an empty object. Please call the build_Minigfa method to construct
    Examples:
            >>> myGfa = Minigfa("pangenome.gfa")
            >>> myGfa2 = Minigfa()
            >>> myGfa2.build_Minigfa("pangenome.gfa")

    """

    S_line_factory = dict[str, _Segment]
    L_line_factory = list[_Link]

    def __init__(self, file_path: Optional[str] = None) -> None:
        self.S_line = self.S_line_factory()
        self.L_line = self.L_line_factory()
        if file_path != None:
            self.build_Minigfa(file_path)

    def build_Minigfa(self, file_path: str) -> None:
        try:
            with open(file_path, "r") as file:
                for lineno, line in enumerate(file, start=1):
                    temp = line.strip()
                    if temp.startswith("S"):
                        easy_line = line.strip().split("\t")
                        self.S_line[easy_line[1]] = _Segment(temp)
                    elif temp.startswith("L"):
                        self.L_line.append(_Link(temp))
                    else:
                        continue
        except FileNotFoundError:
            print(f"Error: File '{file_path}' not found.")
            exit(1)

    def get_linear_reference(self) -> str:
        return self.S_line["s1"].source_sample

    def get_all_segID(self) -> Generator[str, Any, None]:
        """Generate all segment IDs

        Yields:
            str: return a segment ID
        """
        for segID in self.S_line.keys():
            yield segID

    def get_seq(self, segID: str) -> str:
        return self.S_line[segID].seq

    def get_source_sample(self, segID: str) -> str:
        return self.S_line[segID].source_sample

    def get_SRank(self, segID: str) -> int:
        return self.S_line[segID].SRank

    def get_all_Link(self) -> Generator[tuple[str, str, str, str, int], Any, None]:
        """Generate all segment meesages

        Yields:

            tuple[str, str, str, str, int]: Returns a five-tuple, fromID, fromOrient, toID,toOrient, SRank in order
        """
        for link in self.L_line:
            yield link.fromID, link.fromOrient, link.toID, link.toOrient, link.SRank


class _bedLine:
    def __init__(self, bed_line: str) -> None:
        easy_line = bed_line.strip().split("\t")
        self.chr: str = easy_line[0].strip().split("#")[2]
        self.segs_num: int = int(easy_line[3])
        self.possible_paths_num: int = int(easy_line[4])
        self.list_of_segments: list[str] = easy_line[11].strip().split(",")
        self.is_inverved: bool = True if int(easy_line[5]) == 1 else False


class Minibed:
    """
    Construct a composite data storing BED file information.
       This class is iterable.Each iteration yields a five-tuple, chr_num, is_invered, segs_num, possibal_path_num,list_of_segments in order from one line in BED file.

       The path of the BED file that is preferably passed in when constructing the object.If you don't do this, you will just get an empty object. Please call the build_Minigfa method to construct.
    """

    def __init__(self, file_path: Optional[str] = None) -> None:
        self.filePath: Optional[str] = None
        if file_path is not None:
            self.filePath = file_path

    def build_Minibed(self, file_path: str) -> None:
        self.filePath = file_path

    def __iter__(self) -> Generator[tuple[str, bool, int, int, list[str]], Any, None]:
        if self.filePath is None:
            raise TypeError("You must give Minibed a file_path")
        try:
            with open(self.filePath, "r") as fileStream:
                for line in fileStream:
                    temp = _bedLine(line)
                    yield temp.chr, temp.is_inverved, temp.segs_num, temp.possible_paths_num, temp.list_of_segments
        except FileNotFoundError:
            print(f"Error: File '{self.filePath}' not found.")
            exit(1)

    @cache
    def get_linear_sources_and_sinks(self) -> tuple[dict[str, str], dict[str, str]]:
        bed_line = dict[str, list[str]]()
        sources = dict[str, str]()
        sinks = dict[str, str]()
        for chr, _, _, _, list_of_Segment in self:
            if chr in bed_line.keys():
                bed_line[chr].extend(list_of_Segment)
            else:
                bed_line[chr] = list_of_Segment
        for chr_and_bedLine in bed_line.items():
            sources[chr_and_bedLine[0]] = chr_and_bedLine[1][0]
            sinks[chr_and_bedLine[0]] = chr_and_bedLine[1][-1]
        return sources, sinks


if __name__ == "__main__":
    # MiniBed = Minibed("/io/wh/hprc_v4.0/CHM13-464.bb.bed")
    # sources, sinks = MiniBed.get_linear_sources_and_sinks()
    # print(sources)
    # print(sinks)
    pass
