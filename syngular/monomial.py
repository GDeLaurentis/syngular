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

    def __new__(cls, *args):
        if len(args) == 2 and isinstance(args[0], (list, tuple)) and isinstance(args[1], (list, tuple)):
            combined_arg = (dict(zip(args[0], args[1])),)
            return super(Monomial, cls).__new__(cls, *combined_arg)
        return super(Monomial, cls).__new__(cls, *args)

    def __init__(self, *args):
        if len(args) == 0:
            args = [(), ]
        if isinstance(args[0], (dict, FrozenMultiset, Monomial)):
            super(Monomial, self).__init__(args[0])
        elif isinstance(args[0], (tuple, list)) and (len(args) == 1 or not isinstance(args[1], (tuple, list))):
            super(Monomial, self).__init__(self.__rstr__('·'.join(args[0])))
        elif isinstance(args[0], str):
            super(Monomial, self).__init__(self.__rstr__(args[0]))
        elif len(args) >= 2 and isinstance(args[0], (list, tuple)) and isinstance(args[1], (list, tuple)):
            super(Monomial, self).__init__(dict(zip(args[0], args[1])))
        else:
            raise NotImplementedError(
                f"""Monomial initialization not understood, received:\n
                args: {args}\n
                args types: {list(map(type, args))}
                """)

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
        return list(self.keys())

    def subs(self, values_dict):
        return functools.reduce(operator.mul, [values_dict[key] ** val for key, val in self.items()], 1)

    def __mul__(self, other):
        assert isinstance(other, Monomial)
        return super(Monomial, self).__add__(other)

    def __truediv__(self, other):
        if isinstance(other, Monomial) and other <= self:
            return Monomial(FrozenMultiset(self) - FrozenMultiset(other))
        else:
            return NotImplemented
        # raise Exception("Monomial division not implement. Do you mean this to be a Rational Function?")

    def __add__(self, other):
        raise Exception("Monomial addition not implement. Do you mean this to be a Polynomial?")

    def __sub__(self, other):
        raise Exception("Monomial subtraction not implement. Do you mean this to be a Polynomial?")
