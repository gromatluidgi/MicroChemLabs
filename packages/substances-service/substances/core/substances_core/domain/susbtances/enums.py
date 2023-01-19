from enum import Enum


class ReactivityMode(str, Enum):
    """
    Nearly all chemical reactions, whether organic or inorganic,
    proceed because atoms or groups of atoms having a positive charge
    or a partial positive charge interact with atoms or groups of atoms having
    a negative charge or a partial negative charge.

    Thus when a bond in a hydrocarbon is cleaved during a reaction,
    identifying the transient species formed, some of which are charged,
    allows chemists to determine the mechanism and predict the products of a reaction.
    https://chem.libretexts.org/Bookshelves/General_Chemistry/Book%3A_General_Chemistry%3A_Principles_Patterns_and_Applications_(Averill)/23%3A_Organic_Compounds/23.04%3A_Reactivity_of_Organic_Molecules
    """

    ELECTROPHILE = "electron_deficient"
    NUCLEOPHILE = "electron_rich"
    NONE = "none"


class ReviewStatus(str, Enum):
    """ReviewStatus"""

    PENDING = "pending"
    VALID = "valid"
    DECLINED = "declined"
