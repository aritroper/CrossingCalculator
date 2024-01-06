from enum import Enum, auto
from itertools import product
import copy

class Gender(Enum):
    MALE = auto()
    FEMALE = auto()

    def __str__(self):
        if self == Gender.MALE:
            return "(M)"
        elif self == Gender.FEMALE:
            return "(F)"
        else:
            raise ValueError("Gender must be M or F")

class Phenotype(Enum):
    WHITE_EYES = auto()
    ORANGE_EYES = auto()
    RED_EYES = auto()
    STAR_EYES = auto()
    CURLY_WINGS = auto()
    SHORT_BRISTLES = auto()
    BIG_HALTERES = auto()
    SHOULDER_HAIR = auto()
    MALE = auto()
    RECOMBINATED = auto() # Not really a phenotype but can be screened for with PCR
    NONE = auto()

    def __str__(self):
        return ' '.join(word.capitalize() for word in self.name.split('_'))

    # Implement the less than method for sorting
    def __lt__(self, other):
        if isinstance(other, Phenotype):
            return self.value < other.value
        return NotImplemented

class Gene:
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

    def __eq__(self, other):
        if isinstance(other, Gene):
            return self.name == other.name
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, Gene):
            return self.name < other.name
        return NotImplemented

    def __hash__(self):
        return hash(self.name)

class Marker(Gene):
    def __init__(self, name, phenotype):
        super().__init__(name, phenotype, True, True, True)

class Balancer(Gene):
    def __init__(self, name, phenotype):
        super().__init__(name, phenotype, False, True, True)

class Transgene(Gene):
    def __init__(self, name):
        super().__init__(name, Phenotype.NONE, True, False, False)

    def get_phenotype(self, has_white_eyes):
        if has_white_eyes:
            self.visible = True
            return Phenotype.ORANGE_EYES
        else:
            self.visible = False
            return Phenotype.NONE

class Mutation(Gene):
    def __init__(self, name, phenotype):
        super().__init__(name, phenotype, True, False, True)

class GeneType():

    Y = Gene("Y", Phenotype.MALE, False, True, True)
    PLUS = Gene("PLUS", Phenotype.NONE, True, False, False)

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

    # Transgenes
    XGAL4 = Transgene("XGAL4")
    DNAD = Transgene("DN-AD")
    SOS = Transgene("SOS")
    DNDB = Transgene("DN-DB")
    ORB = Transgene("ORB")
    UAS = Transgene("UAS-SPARC")

    # Mutations
    WMINUS = Mutation("W-", Phenotype.WHITE_EYES)
    WPLUS = Mutation("W+", Phenotype.RED_EYES)

class Chromatid():
    def __init__(self, genes):
        if len(genes) > 2:
            raise ValueError("Chromatid should have max two genes on it")
        self.genes = genes

    '''
    Called when sister chromatid is not specified in the chromosome
    '''
    def get_implied_sister_chromatid(self):
        if not self.genes:
            raise ValueError("Chromatid must have at least one gene to get implied sister chromatid")

        return Chromatid([GeneType.PLUS if gene.homozygous_lethal else gene for gene in self.genes])

    def get_phenotype(self, has_white_eyes):
        return {gene.get_phenotype(has_white_eyes) for gene in self.genes}

    def num_genes(self):
        return len(self.genes)

    def get_visual_genes(self):
        return {gene for gene in self.genes if gene.visible}

    '''
    Given a list of genes, returns the first index
    to match one of these genes
    '''
    def get_first_idx_of(self, genes):
        for i in range(len(genes)):
            for gene in genes:
                if genes[i] == gene:
                    return i

    def __getitem__(self, index):
        """Get the gene at the specified index."""
        if 0 <= index < len(self.genes):
            return self.genes[index]
        else:
            raise IndexError("Index out of range")

    def __str__(self):
        return "[" + ' '.join(gene.name for gene in self.genes) + "]"

    def __eq__(self, other):
        return sorted(self.genes) == sorted(other.genes)

    def __hash__(self):
        # Sort the genes and concatenate their names to create a unique string
        sorted_gene_names = ''.join(sorted(gene.name for gene in self.genes))
        return hash(sorted_gene_names)

class Chromosome():
    def __init__(self, chromatid1, chromatid2=None):

        # Casting magic
        if isinstance(chromatid1, Gene):
            chromatid1 = Chromatid([chromatid1])

        if isinstance(chromatid2, Gene):
            chromatid2 = Chromatid([chromatid2])

        if isinstance(chromatid1, Chromatid) and chromatid2 is None:
            self.chromatid1 = chromatid1
            self.chromatid2 = chromatid1.get_implied_sister_chromatid()
        elif isinstance(chromatid2, Chromatid) and chromatid1 is None:
            self.chromatid1 = chromatid2.get_implied_sister_chromatid()
            self.chromatid2 = chromatid2
        else:
            self.chromatid1 = chromatid1
            self.chromatid2 = chromatid2

        for gene in self.chromatid1.genes:
            if gene.homozygous_lethal and gene in self.chromatid2.genes:
                raise ValueError("Cannot construct this chromosome")

    def get_phenotype(self, has_white_eyes):
        return sorted(chromatid1.get_phenotype(has_white_eyes) + chromatid2.get_phenotype(has_white_eyes))

    def can_undergo_recombination(self):
        # Limit number of recombinations to 1 (for now...)
        # TODO W+ W-
        if self.chromatid1.num_genes() == 1 and self.chromatid2.num_genes() == 1:
            return self.chromatid1[0].recombination and self.chromatid2[0].recombination and (self.chromatid1[0] != self.chromatid2[0])
        else:
            return False

    '''
    Returns recombinated copy of self
    '''
    def recombinate(self):
        # TODO double check with Janina
        # To make this more efficient the chromosome does not need to be deep copied
        # since the line is deep copied when recombinating line
        if self.can_undergo_recombination():
            recombinated = copy.deepcopy(self)
            recombinated.chromatid1.genes.append(recombinated.chromatid2[0])
            recombinated.chromatid2.genes[0] = GeneType.PLUS
            return recombinated
        else:
            raise ValueError("Chromosome cannot be recombinated")

    '''
    Returns True if this chromosome has undergone at least 1 recombination
    '''
    def has_recombinated(self):
        return len(self.chromatid1.genes) > 1 or len(self.chromatid2.genes) > 1

    # IMPORTANT: Should only be called on the first chromosome
    def has_white_eyes(self):
        has_wminus_in_chromatid1 = GeneType.WMINUS in self.chromatid1.genes
        has_wminus_in_chromatid2 = GeneType.WMINUS in self.chromatid2.genes

        has_y_in_chromatid1 = GeneType.Y in self.chromatid1.genes
        has_y_in_chromatid2 = GeneType.Y in self.chromatid2.genes

        return (has_wminus_in_chromatid1 and (has_wminus_in_chromatid2 or has_y_in_chromatid2)) or \
               (has_wminus_in_chromatid2 and (has_wminus_in_chromatid1 or has_y_in_chromatid1))

    # IMPORTANT: Should only be called on the first chromosome
    def is_male(self):
        return (GeneType.Y in self.chromatid1.genes) or (GeneType.Y in self.chromatid2.genes)

    def get_visual_genes(self):
        return self.chromatid1.get_visual_genes() | self.chromatid2.get_visual_genes()

    def get_phenotype(self):
        return self.chromatid1.get_visual_genes() | self.chromatid2.get_visual_genes()

    def generate_haploids(self):
        return [self.chromatid1.genes, self.chromatid2.genes]

    def __str__(self):
        if (self.chromatid1 == self.chromatid2):
            return f"{str(self.chromatid1)}"
        else:
            return f"{str(self.chromatid1)}/{str(self.chromatid2)}"

    def __eq__(self, other):
        return ((self.chromatid1, self.chromatid2) == (other.chromatid1, other.chromatid2)) or ((self.chromatid2, self.chromatid1) == (other.chromatid1, other.chromatid2))

    def __hash__(self):
        # TODO Potentially review
        return hash(self.chromatid1) + hash(self.chromatid2)

class Line:
    def __init__(self, chromosome1=GeneType.PLUS, chromosome2=GeneType.PLUS, chromosome3=GeneType.PLUS, chromosome4=GeneType.PLUS, id=-1):
        genotype = [chromosome1, chromosome2, chromosome3, chromosome4]

        for i in range(len(genotype)):
            # Cast genes to chromosomes
            if isinstance(genotype[i], Gene):
                genotype[i] = Chromosome(Chromatid([genotype[i]]))

        self.genotype = genotype
        self.has_white_eyes = self[0].has_white_eyes()
        self.id = id
        self.has_recombinated = any(chromosome.has_recombinated() for chromosome in self.genotype)
        
        if (self[0].is_male()):
            self.gender = Gender.MALE
        else:
            self.gender = Gender.FEMALE

    def get_genotype(self):
        """Return a list of genes in the starter line."""
        return self.genotype

    def get_visual_genes(self):
        return set().union(*(chromosome.get_visual_genes() for chromosome in self.genotype))

    def can_undergo_recombination(self):
        return (self.gender) == Gender.FEMALE and (not self.has_recombinated) and (any(chromosome.can_undergo_recombination() for chromosome in self.genotype))

    '''
    Returns recombinated copy of self
    '''
    def recombinate(self):
        if self.can_undergo_recombination():
            recombinated = copy.deepcopy(self)
            for i in range(4):
                if recombinated.genotype[i].can_undergo_recombination():
                    recombinated.genotype[i] = recombinated.genotype[i].recombinate()
            recombinated.has_recombinated = True
            return recombinated
        else:
             raise ValueError("Line cannot be recombinated")

    def get_phenotype(self):
        # TODO ask janina
        phenotypes = set().union(*(chromosome.get_phenotype() for chromosome in self.genotype))
        if self.has_recombinated:
            # Lines that have recombinated can be selected for with PCR screening even if same phenotype
            phenotypes.add(Phenotype.RECOMBINATED)

        return phenotypes

    def flip_gender(self):
        opposite_gender_line = copy.deepcopy(self)
        if self.gender == Gender.MALE:
            # Replace Y with duplicate copy of W- or W+
            opposite_gender_line.gender = Gender.FEMALE
            idx_to_replace = opposite_gender_line[0].chromatid2.get_first_idx_of([GeneType.Y])
            opposite_gender_line[0].chromatid2.genes[idx_to_replace] = opposite_gender_line[0].chromatid1.genes[idx_to_replace]
        else:   
            # Replace W- or W+ with Y
            opposite_gender_line.gender = Gender.MALE
            idx_to_replace = opposite_gender_line[0].chromatid2.get_first_idx_of([GeneType.WMINUS, GeneType.WPLUS])
            opposite_gender_line[0].chromatid2.genes[idx_to_replace] = GeneType.Y
        return opposite_gender_line

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
        chrom_strs = []
        for chromosome in self.genotype:
            if isinstance(chromosome, Chromosome):
                chrom_strs.append(str(chromosome))
            else:
                chrom_strs.append("??")  # Fallback for unexpected types

        return ", ".join(chrom_strs)


def cross_chromosome(chromosome1, chromosome2):
        if not isinstance(chromosome1, Chromosome) or not isinstance(chromosome2, Chromosome):
            raise ValueError("cross_chromosome Error: both chromosomes should be represented as chromosomes")

        crossed = set()  # Initialize an empty set
        chromatids1 = chromosome1.generate_haploids()
        chromatids2 = chromosome2.generate_haploids()

        # Perform initial crossing
        for chromatid1 in chromatids1:
            for chromatid2 in chromatids2:
                try:
                    c = Chromosome(Chromatid(chromatid1), Chromatid(chromatid2))
                    crossed.add(c)
                except ValueError:
                    continue

        return crossed

def _cross_lines(male_line, female_line):
    # Must cross Male with Female
    if male_line.gender != Gender.MALE or female_line.gender != Gender.FEMALE:
        raise ValueError("Must cross male line with female line")

    all_chromosome_gene_combinations = []

    for i in range(4):
        crossed_genes = cross_chromosome(female_line[i], male_line[i])
        all_chromosome_gene_combinations.append(crossed_genes)

    return [Line(*combination) for combination in product(*all_chromosome_gene_combinations)]

def cross_lines(line1, line2):
    # Must cross Male with Female
    if line1.gender == line2.gender:
        return []

    male_line = line1 if line1.gender == Gender.MALE else line2
    female_line = line1 if line1.gender == Gender.FEMALE else line2

    lines_no_recombination = _cross_lines(male_line, female_line)

    if (female_line.can_undergo_recombination()):
        lines_with_recombination = _cross_lines(male_line, female_line.recombinate())
        lines = lines_no_recombination + lines_with_recombination
    else:
        lines = lines_no_recombination

    # Handling non-hashable phenotypes
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

def compute_crosses(starter_lines, target):
    seen_lines = set(starter_lines)
    queue = [(line, []) for line in starter_lines]

    while queue:
        new_queue = [(line, []) for line in starter_lines]
        added_any_new_line = False  # Flag to check if new lines are added

        for i, (line1, path1) in enumerate(queue):
            for j in range(i + 1, len(queue)):  # Start from i + 1
                line2, path2 = queue[j]
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

    return None

def compute_pick_for(line1, line2, progeny_line):
    genes1 = line1.get_visual_genes()
    genes2 = line2.get_visual_genes()

    all_genes = genes1 | genes2

    pick_for = set()
    pick_against = set()

    progeny_genes = progeny_line.get_visual_genes()

    # Compute pick-for, pick-against
    for gene in all_genes:
        if gene in progeny_genes:
            pick_for.add(gene)
        else:
            pick_against.add(gene)

    # Print pick-for, pick-against
    for g in pick_for:
        print("+" + str(g))

    for g in pick_against:
        print("-" + str(g))

def init_starter_lines(starter_lines):
    LINE_ID = 1
    all_lines = []

    # Add both Male and Female versions of each starting line
    for line in starter_lines:
        line_gender_flipped = line.flip_gender()

        line.id = LINE_ID
        line_gender_flipped.id = LINE_ID

        all_lines.append(line)
        all_lines.append(line_gender_flipped)

        LINE_ID += 1

    return all_lines, LINE_ID

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
        print("(" + str(lines[0].id) + ") " + str(lines[0]) + " | " + str(lines[0].gender))
        print(" X")
        print("(" + str(lines[1].id) + ") " + str(lines[1]) + " | " + str(lines[1].gender))
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
lineA = Line(GeneType.WMINUS, GeneType.PLUS, GeneType.PLUS, GeneType.ORB)
lineB = Line(GeneType.WMINUS, GeneType.PLUS, GeneType.DNAD, GeneType.PLUS)
lineC = Line(GeneType.WMINUS, GeneType.XGAL4, GeneType.PLUS, GeneType.PLUS)

# Target
target = Line(GeneType.WMINUS, GeneType.XGAL4, GeneType.DNAD, GeneType.ORB)

main([lineA, lineB, lineC], target)
