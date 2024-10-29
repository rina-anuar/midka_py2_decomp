import re

class DerivativeOfPowerFunctions:
    def __init__(self, equation):
        self.equation = equation
        self.dydx = ""
        self.is_polynomial = False
        self.is_trigonometric = False

    def derivative_terms(self, term):
        coeff = None
        pow = None

        # Поиск коэффициента и степени
        if re.search(r"\*x", term):
            coeff = int(re.findall(r"([+-]?\d*)\*x", term)[-1] or '1')  # если нет коэффициента, предполагаем 1
            pow = int(re.findall(r"\^(\d+)", term)[-1]) if re.search(r"\^", term) else 1  # степень по умолчанию 1

            # Вычисление производной
            if pow is not None:
                if pow != 1:
                    return f"{coeff * pow:+}*x^{pow - 1}" if pow != 2 else f"{coeff * pow:+}*x"
                else:
                    return f"{coeff:+}"
        return "0"  # если не найдено x, производная равна 0

    def diff_polynomial(self):
        terms = re.findall(r"([+-]?\d*\*x\^?\d*)", self.equation)
        self.dydx = ""

        for term in terms:
            self.dydx += self.derivative_terms(term)

        return self.dydx.strip("+")

    def diff(self):
        if len(re.findall(r"^_?\d+$", self.equation)) != 0:
            return "0"
        elif len(re.findall(r"^([+-]?\d*\*x\^?\d*)", self.equation)) != 0:
            self.is_polynomial = True
            return self.diff_polynomial()
        elif len(re.findall(r"-?sin\(x\)|-?cos\(x\)|-?e\^x|-?ln\(x\)", self.equation)) != 0:
            self.is_trigonometric = True
            print("Trigonometric differentiation:", self.equation)
            # Здесь нужно добавить логику для тригонометрических функций
            return self.dydx  # вернем текущую производную
        return "0"

# Пример использования
equation = "3*sin(x)"
derivative_calculator = DerivativeOfPowerFunctions(equation)
print(derivative_calculator.diff())