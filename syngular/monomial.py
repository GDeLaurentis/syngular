import functools
import operator
import re
import syngular

from multiset import FrozenMultiset

unicode_powers_dict = {"0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
                       "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹"}


def unicode_powers(string):
    """Convert ^ followed by digits into superscripts."""
    return re.sub(r"\^(\d+)", lambda m: "".join(unicode_powers_dict[d] for d in m.group(1)), string)


def non_unicode_powers(string):
    """Convert superscript digits back into ^ followed by normal digits."""
    return re.sub(r"[⁰¹²³⁴⁵⁶⁷⁸⁹]+",
                  lambda m: "^" + "".join(str(list(unicode_powers_dict.values()).index(c)) for c in m.group(0)),
                  string)


class Monomial(FrozenMultiset):
    """A FrozenMultiset representation of a Monomial. Positive integer multiplicities represent powers."""

    def __init__(self, letters_and_powers={}):
        if isinstance(letters_and_powers, (dict, tuple, FrozenMultiset, Monomial)):
            super(Monomial, self).__init__(letters_and_powers)
        elif isinstance(letters_and_powers, str):
            super(Monomial, self).__init__(self.__rstr__(letters_and_powers))
        else:
            print("entry:", repr(letters_and_powers))
            print("type:", type(letters_and_powers))
            raise NotImplementedError

    @staticmethod
    def __rstr__(string):
        string = non_unicode_powers(string).replace("**", "^")
        while len(string) > 0 and string[0] == ' ':
            string = string[1:]
        if string == '':
            return dict()
        splitted_string = [entry for entry in re.split(r"(?<![\+\-\(])(?<!tr5)([⟨\[]|(?<![a-zA-Z])[\(a-zA-ZΔΩΠ])", string) if entry != '']
        splitted_string = [splitted_string[i] + splitted_string[i + 1] if i + 1 < len(splitted_string) else splitted_string[i]
                           for i in range(len(splitted_string))[::2]]
        # sqeuentially remerge strings until parenthesis are (minimally) balanced
        splitted_string_partially_remerged = [splitted_string[0]]
        for entry in splitted_string[1:]:
            if splitted_string_partially_remerged[-1].count("(") != splitted_string_partially_remerged[-1].count(")"):
                splitted_string_partially_remerged[-1] += entry
            else:
                splitted_string_partially_remerged += [entry]
        splitted_string_partially_remerged = [inv_and_exp[:-1] if inv_and_exp[-1] in ['*', '·', ' '] else inv_and_exp
                                              for inv_and_exp in splitted_string_partially_remerged]
        # print(splitted_string_partially_remerged)
        invs, exps = list(map(list, zip(*[re.split(r"\^(?=\d+$)", entry) if re.findall(r"\^(?=\d+$)", entry) != []
                                        else [entry, '1'] for entry in splitted_string_partially_remerged])))
        invs = [inv.replace(" ", "") for inv in invs]
        # print(list(zip(invs, list(map(int, exps)))))
        return dict(zip(invs, list(map(int, exps))))

    def __repr__(self):
        return f"Monomial(\"{str(self)}\")"

    def __str__(self):
        string = re.sub(r"\^1(?![\.\d])", "", f"{syngular.CDOTCHAR}".join([f"{key}^{val}" for key, val in self.items()]))
        if syngular.UNICODEPOWERS:
            string = unicode_powers(string)
        if not syngular.FORCECDOTS:
            string = re.sub(rf'(?<![a-zA-Z]){syngular.CDOTCHAR}|{syngular.CDOTCHAR}(?![a-zA-Z])', '', string)
        return string

    def tolist(self):
        return [entry for key, val in self.items() for entry in [key, ] * val]

    def subs(self, values_dict):
        return functools.reduce(operator.mul, [values_dict[key] ** val for key, val in self.items()], 1)

    def __mul__(self, other):
        assert isinstance(other, Monomial)
        return super(Monomial, self).__add__(other)

    def __truediv__(self, other):
        raise Exception("Monomial division not implement. Do you mean this to be a Rational Function?")

    def __add__(self, other):
        raise Exception("Monomial addition not implement. Do you mean this to be a Polynomial?")

    def __sub__(self, other):
        raise Exception("Monomial subtraction not implement. Do you mean this to be a Polynomial?")
