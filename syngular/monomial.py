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


def preserve_class_binary_op(func):
    @functools.wraps(func)
    def wrapper(self, other):
        result = func(self, other)
        if result is NotImplemented:
            return NotImplemented
        cls = self.__class__ if issubclass(self.__class__, other.__class__) else other.__class__
        return cls(result)
    return wrapper


class Monomial(FrozenMultiset):
    """A FrozenMultiset representation of a Monomial. Positive integer multiplicities represent powers."""

    def __new__(cls, *_):
        return super(Monomial, cls).__new__(cls, )

    def __init__(self, *args):
        if len(args) == 0:
            args = [(), ]
        if isinstance(args[0], (dict, FrozenMultiset, Monomial)):
            assert all([isinstance(exp, int) or (isinstance(exp, float) and exp.is_integer()) for exp in args[0].values()])
            data = {key: int(val) for key, val in args[0].items()}
        elif len(args) == 1 and isinstance(args[0], (tuple, list)) and all([isinstance(entry, (list, tuple)) for entry in args[0]]):
            data = [inv for inv, exp in args[0] for _ in range(int(exp))]
        elif len(args) == 1 and isinstance(args[0], (tuple, list)) and all([isinstance(entry, str) for entry in args[0]]):
            data = self.__rstr__('·'.join(args[0]))
        elif len(args) == 2 and isinstance(args[0], (list, tuple)) and isinstance(args[1], (list, tuple)):
            assert all([isinstance(exp, int) or (isinstance(exp, float) and exp.is_integer()) for exp in args[1]])
            data = [inv for inv, exp in zip(args[0], map(int, args[1])) for _ in range(exp)]
        elif isinstance(args[0], str):
            data = self.__rstr__(args[0])
        else:
            raise NotImplementedError(f"Monomial initialization not understood, received:\nargs: {args}\nargs types: {list(map(type, args))}")
        for pattern in syngular.NORMALIZE_POWERS_PATTERNS:
            pattern_complete = re.compile(rf"{pattern.pattern}(\^\d+)?")
            pattern_full_match = re.compile(rf"^{pattern.pattern}$")
            data = [(pattern_complete.sub(
                lambda match: f"{match.group(1)}^{match.group(2)}" if match.group(3) is None else
                              f"{match.group(1)}^{int(match.group(2)) * int(match.group(3)[1:])}",
                              inv), exp) if pattern_full_match.findall(inv) == [] else
                    (pattern_full_match.findall(inv)[0][0], int(pattern_full_match.findall(inv)[0][1]) * exp)
                    for inv, exp in FrozenMultiset(data).items()]
            data = [inv for inv, exp in data for _ in range(exp)]
        super(Monomial, self).__init__(data)

    @staticmethod
    def __rstr__(string):
        # print(string)
        string = non_unicode_powers(string).replace("**", "^")
        string = " ".join(string.split())
        if string == '':
            return dict()
        pattern = r"(?<![\+\-\(])(?<!tr5)(?<!tr)(?=[⟨\[\(]|(?<![a-zA-ZΔΩΠ_])[a-zA-ZΔΩΠ])"
        splitted_string = [entry for entry in re.split(pattern, string) if entry]
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
        invs = [' '.join(inv.split()) for inv in invs]
        # print(list(zip(invs, list(map(int, exps)))))
        return FrozenMultiset(inv for inv, exp in zip(invs, exps) for _ in range(int(exp)))

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

    def __call__(self, other):
        if callable(other):
            return self.subs({key: other(key) for key in self.keys()})
        return self.subs(other)

    @preserve_class_binary_op
    def __mul__(self, other):
        if isinstance(other, Monomial):
            return FrozenMultiset(self) + FrozenMultiset(other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return other * self

    @preserve_class_binary_op
    def __truediv__(self, other):
        if isinstance(other, Monomial) and other.issubset(self):
            return FrozenMultiset(self) - FrozenMultiset(other)
        else:
            return NotImplemented
        # raise Exception("Monomial division not implement. Do you mean this to be a Rational Function?")

    def __add__(self, other):
        raise Exception("Monomial addition not implement. Do you mean this to be a Polynomial?")

    def __sub__(self, other):
        raise Exception("Monomial subtraction not implement. Do you mean this to be a Polynomial?")

    @property
    def invs(self):
        return [inv for inv, _ in self.items()]

    @property
    def exps(self):
        return [exp for _, exp in self.items()]

    @property
    def variables(self):
        return set(self.invs)

    def __lt__(self, other):
        # Graded reverse lexicographic ordering
        self_deg = sum(self.values())
        other_deg = sum(other.values())
        if self_deg != other_deg:
            return self_deg > other_deg  # higher total degree comes first

        self_keys = sorted(self.keys(), reverse=True)
        other_keys = sorted(other.keys(), reverse=True)
        for k1, k2 in zip(self_keys, other_keys):
            if k1 != k2:
                return k1 < k2
            if self[k1] != other[k2]:
                return self[k1] < other[k2]
        return len(self) < len(other)

    def __le__(self, other):
        return self == other or self < other

    def __gt__(self, other):
        return not self <= other

    def __ge__(self, other):
        return not self < other
