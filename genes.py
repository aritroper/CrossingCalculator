from enum import Enum, auto
from itertools import product

class Gender(Enum):
    MALE = auto()
    FEMALE = auto()
    BOTH = auto()

    def __str__(self):
        if self == Gender.MALE:
            return "(M)"
        elif self == Gender.FEMALE:
            return "(F)"
        else:
            return "(M/F)"

class Phenotype(Enum):
    WHITE_EYES = auto()
    ORANGE_EYES = auto()
    RED_EYES = auto()
    STAR_EYES = auto()
    CURLY_WINGS = auto()
    SHORT_BRISTLES = auto()
    BIG_HALTERES = auto()
    SHOULDER_HAIR = auto()
    NONE = auto()

    def __str__(self):
        return ' '.join(word.capitalize() for word in self.name.split('_'))

    # Implement the less than method for sorting
    def __lt__(self, other):
        if isinstance(other, Phenotype):
            return self.value < other.value
        return NotImplemented

class Allele:
    def __init__(self, name, phenotype, recombination, homozygous_lethal, visible):
        self.name = name
        self.phenotype = phenotype
        self.recombination = recombination
        self.homozygous_lethal = homozygous_lethal
        self.visible = visible

    def __str__(self):
        return self.name

    def get_phenotype(self, has_white_eyes):
        return self.phenotype

class Marker(Allele):
    def __init__(self, name, phenotype):
        super().__init__(name, phenotype, True, True, True)

class Balancer(Allele):
    def __init__(self, name, phenotype):
        super().__init__(name, phenotype, False, True, True)

class Transgene(Allele):
    def __init__(self, name, mwc):
        super().__init__(name, Phenotype.NONE, True, False, False)
        self.mwc = mwc

    def get_phenotype(self, has_white_eyes):
        if has_white_eyes and mwc:
            self.visible = True
            return Phenotype.ORANGE_EYES
        else:
            self.visible = False
            return Phenotype.NONE

class Mutation(Allele):
    def __init__(self, name, phenotype):
        super().__init__(name, phenotype, True, False, True)

class AlleleType():

    Y = Allele("Y", "Male", False, True, True)
    PLUS = Allele("PLUS", Phenotype.NONE, True, False, False)

    # Markers
    SM6 = Marker("SM6", Phenotype.CURLY_WINGS)
    S = Marker("S", Phenotype.STAR_EYES)
    UBX = Marker("UBX", Phenotype.BIG_HALTERES)
    TB = Marker("TB", Phenotype.SHOULDER_HAIR)

    # Balancers
    CyO = Balancer("CyO", Phenotype.CURLY_WINGS)
    TM2 = Balancer("TM2", Phenotype.BIG_HALTERES)
    TM3 = Balancer("TM3", Phenotype.SHORT_BRISTLES)
    TM6b = Balancer("Tm6b", Phenotype.SHOULDER_HAIR)

    # Transgenes (we should add mwc as a parameter)
    XGAL4 = Transgene("XGAL4", True)
    DNAD = Transgene("DN-AD", True)
    SOS = Transgene("SOS", True)
    DNDB = Transgene("DN-DB", True)
    ORB = Transgene("ORB", True)
    UAS = Transgene("UAS-SPARC", True)

    # Mutations
    WMINUS = Mutation("W-", Phenotype.WHITE_EYES)
    WPLUS = Mutation("W+", Phenotype.RED_EYES)

class Gene():
    def __init__(self, allele1, allele2=None):
        if allele2 is None:
            # Only one gene provided, set gene2 to a default value
            if allele1.homozygous_lethal:
                allele2 = AlleleType.PLUS
            else:
                allele2 = allele1
        if (allele1 == allele2) and (allele1.homozygous_lethal or allele2.homozygous_lethal):
            raise ValueError("Gene Error: allele1 and allele2 must be different")
        if not isinstance(allele1, Allele) or not isinstance(allele2, Allele):
            raise ValueError("Gene Error: allele1 and allele2 must be instances of allele")

        self.allele1 = allele1
        self.allele2 = allele2

    def get_allele_not_transgene(self):
        if isinstance(self.allele1, Transgene) and isinstance(self.allele2, Transgene):
            return None
        elif isinstance(self.allele1, Transgene):
            return self.allele2
        else:
            return self.allele1

    def get_phenotype(self, has_white_eyes):
        return sorted([self.allele1.get_phenotype(has_white_eyes), self.allele2.get_phenotype(has_white_eyes)])

    def __str__(self):
        if (self.allele1 == self.allele2):
            return f"{self.allele1.name}"
        elif (self.allele1 == AlleleType.PLUS and self.allele2 != AlleleType.PLUS):
            return f"{self.allele2.name}/{self.allele1.name}"
        else:
            return f"{self.allele1.name}/{self.allele2.name}"

    def __eq__(self, other):
        return ((self.allele1, self.allele2) == (other.allele1, other.allele2)) or ((self.allele2, self.allele1) == (other.allele1, other.allele2))

    def __hash__(self):
        return hash(self.allele1) + hash(self.allele2)

class Line:
    def __init__(self, allele1=AlleleType.PLUS, allele2=AlleleType.PLUS, allele3=AlleleType.PLUS, allele4=AlleleType.PLUS, id=-1):
        genotype = [allele1, allele2, allele3, allele4]

        for i in range(len(genotype)):
            # Check if the gene is an instance of Allele and has homozygous_lethal=True
            if isinstance(genotype[i], Allele):
                genotype[i] = Gene(genotype[i])

        self.genotype = genotype
        self.has_white_eyes = (self[0] == Gene(AlleleType.WMINUS, AlleleType.WMINUS)) or (self[0] == Gene(AlleleType.WMINUS, AlleleType.Y))
        self.id = id
        self.gender = Gender.BOTH

    def get_genotype(self):
        """Return a list of genes in the starter line."""
        return self.genotype

    def get_visual_alleles(self):
        alleles = set()
        for i in range(len(self.genotype)):
            if not isinstance(self[i], Gene):
                raise ValueError("All genes in genotype must be of type gene")
            else:
                if (self[i].allele1.visible):
                    alleles.add(self[i].allele1)
                if (self[i].allele2.visible):
                    alleles.add(self[i].allele2)
        return alleles

    def get_phenotype(self):
        phenotype = []
        for gene in self.genotype:
            phenotype.append(gene.get_phenotype(self.has_white_eyes))
        return phenotype

    def __getitem__(self, index):
        """Get the gene at the specified index."""
        if 0 <= index < len(self.genotype):
            return self.genotype[index]
        else:
            raise IndexError("Index out of range")

    def __eq__(self, other):
        return (self[0] == other[0] and self[1] == other[1] and self[2] == other[2] and self[3] == other[3])

    def __hash__(self):
        return (hash(self[0])) + (hash(self[1])) + (hash(self[2])) + (hash(self[3]))

    def __str__(self):
        """String representation of the starter line."""
        gene_strs = []
        for gene in self.genotype:
            if isinstance(gene, Gene):
                gene_strs.append(str(gene))
            else:
                gene_strs.append("Unknown")  # Fallback for unexpected types

        return ", ".join(gene_strs) + " | " + str(self.gender)


def cross_chromosome(gene1, gene2):
        if not isinstance(gene1, Gene) or not isinstance(gene2, Gene):
            raise ValueError("cross_chromosome Error: all genes should be represented as two alleles")

        crossed = set()  # Initialize an empty set
        alleles1 = [gene1.allele1, gene1.allele2]
        alleles2 = [gene2.allele1, gene2.allele2]

        # Perform initial crossing
        for a1 in alleles1:
            for a2 in alleles2:
                try:
                    crossed.add(Gene(a1, a2))  # Use add() to add elements to the set
                except ValueError:
                    continue

        return crossed

def cross_lines(line1, line2):
    all_chromosome_gene_combinations = []

    # crosses each chromosome, adds to list
    for i in range(4):
        crossed_genes = cross_chromosome(line1[i], line2[i])
        all_chromosome_gene_combinations.append(crossed_genes)

    # adds together chromosomes --> generates all possible lines
    lines = [Line(*combination) for combination in product(*all_chromosome_gene_combinations)]

    # Handling non-hashable phenotypes
    # adds unique chromosomes to list
    phenotype_dict = {}
    for line in lines:
        phenotype = line.get_phenotype()
        phenotype_key = str(phenotype)  # Convert phenotype to a string or other hashable form

        if phenotype_key in phenotype_dict:
            phenotype_dict[phenotype_key] = None  # Mark phenotype as duplicated
        elif phenotype_key not in phenotype_dict:
            phenotype_dict[phenotype_key] = line

    # Filter out lines with duplicated phenotypes
    unique_lines = [line for line in phenotype_dict.values() if line is not None]

    return unique_lines

'''
Returns False if the target contains alleles not present in its starter lines.
'''
def dummy_check(starter_lines, target):
     # Collect all alleles from the starter lines
    starter_alleles = set()
    for line in starter_lines:
        for gene in line.genotype:
            starter_alleles.add(gene.allele1)
            starter_alleles.add(gene.allele2)

    # Check if all alleles in the target line are present in the starter lines
    for gene in target.genotype:
        if gene.allele1 not in starter_alleles or gene.allele2 not in starter_alleles:
            return False  # Target contains an allele not found in starter lines

    return True

def compute_crosses(starter_lines, target):
    seen_lines = set(starter_lines)
    queue = [(line, []) for line in starter_lines]

    if not (dummy_check(starter_lines, target)):
        return None

    while queue:
        new_queue = []
        added_any_new_line = False  # Flag to check if new lines are added

        for i, (line1, path1) in enumerate(queue):
            for j, (line2, path2) in enumerate(queue):
                if i == j:  # Skip crossing a line with itself
                    continue

                crossed_lines = cross_lines(line1, line2)

                for crossed_line in crossed_lines:
                    if crossed_line in seen_lines:
                        continue  # Skip if line has already been seen

                    if crossed_line == target:
                        return path1 + path2 + [(line1, line2, crossed_line)]

                    new_queue.append((crossed_line, path1 + path2 + [(line1, line2, crossed_line)]))
                    seen_lines.add(crossed_line)
                    added_any_new_line = True

        if not added_any_new_line:
            return None  # No new lines added, target line is unattainable

        queue = new_queue

def compute_pick_for(line1, line2, end_line):
    alleles1 = line1.get_visual_alleles()
    alleles2 = line2.get_visual_alleles()

    all_alleles = alleles1.union(alleles2)

    pick_for = set()
    pick_against = set()

    alleles_end = end_line.get_visual_alleles()

    # Compute pick-for, pick-against
    for allele in all_alleles:
        if allele in alleles_end:
            pick_for.add(allele)
        else:
            pick_against.add(allele)

    # Print pick-for, pick-against
    for a in pick_for:
        print("+" + str(a))

    for a in pick_against:
        print("-" + str(a))

def init_starter_lines(starter_lines):
    LINE_ID = 1
    for line in starter_lines:
        line.id = LINE_ID
        LINE_ID += 1
    return (starter_lines, LINE_ID)

def print_crosses(crosses, num_starter_lines):
    if crosses == None:
        print("NO SOLUTION FOUND")
        return

    CROSS = 1
    LINE_ID = num_starter_lines

    print(" ")

    for lines in crosses:

        # Assign these lines ids 
        for i in range(3):
            if (lines[i].id == -1):
                lines[i].id = LINE_ID
                LINE_ID += 1

        print("CROSS " + str(CROSS))
        print("(" + str(lines[0].id) + ") " + str(lines[0]))
        print(" X")
        print("(" + str(lines[1].id) + ") " + str(lines[1]))
        print(" ---------------------")
        print("(" + str(lines[2].id) + ") " + str(lines[2]))
        print(" ")

        compute_pick_for(lines[0], lines[1], lines[2])

        print(" ")

        CROSS += 1

def main(starter_lines, target):
    s_lines = init_starter_lines(starter_lines)
    path = compute_crosses(s_lines[0], target)
    print_crosses(path, s_lines[1])

# Starter lines
lineA = Line(AlleleType.WMINUS, Gene(AlleleType.UAS, AlleleType.PLUS), Gene(AlleleType.PLUS, AlleleType.PLUS), AlleleType.PLUS)
lineB = Line(AlleleType.WPLUS, Gene(AlleleType.DNAD, AlleleType.TM3), Gene(AlleleType.PLUS, AlleleType.PLUS), AlleleType.PLUS)
lineC = Line(AlleleType.WPLUS, Gene(AlleleType.PLUS, AlleleType.TM2), Gene(AlleleType.SOS, AlleleType.TM6b), AlleleType.PLUS)
lineD = Line(AlleleType.WMINUS, Gene(AlleleType.PLUS, AlleleType.PLUS), Gene(AlleleType.DNDB, AlleleType.CyO), AlleleType.PLUS)

# Target
target = Line(AlleleType.WMINUS, Gene(AlleleType.DNAD, AlleleType.UAS), Gene(AlleleType.SOS, AlleleType.DNDB), AlleleType.PLUS)

# Compute
main([lineA, lineB, lineC, lineD], target)



