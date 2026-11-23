import hashlib
from fractions import Fraction
from itertools import combinations, product
import json
from math import comb, prod
from pathlib import Path


VECTORS = {'p': (48, 48, 2, 2), 'q': (50, 46, 3, 1)}
EXPECTED = {
    'p': [0, 509608000, 7874904000, 1560192000, 55296000],
    'q': [0, 518463220, 7873379340, 1566757440, 41400000],
}
DENOMINATOR = 100 ** 5


def require(condition, message):
    if not condition:
        raise ValueError(message)


def enumerate_masses(vector):
    masses = [0] * 5
    for word in product(range(4), repeat=5):
        masses[len(set(word))] += prod(vector[label] for label in word)
    return masses


def inclusion_exclusion(vector):
    masses = [0] * 5
    labels = tuple(range(4))
    for count in range(1, 5):
        for support in combinations(labels, count):
            for subset_size in range(count + 1):
                for subset in combinations(support, subset_size):
                    masses[count] += (-1) ** (count - subset_size) * sum(
                        vector[label] for label in subset
                    ) ** 5
    return masses


def multinomial_masses(vector):
    masses = [0] * 5
    for first in range(6):
        for second in range(6 - first):
            for third in range(6 - first - second):
                fourth = 5 - first - second - third
                counts = (first, second, third, fourth)
                coefficient = comb(5, first) * comb(5 - first, second) * comb(
                    5 - first - second, third
                )
                masses[sum(count > 0 for count in counts)] += coefficient * prod(
                    weight ** count for weight, count in zip(vector, counts)
                )
    return masses


def polynomial(coefficients):
    values = list(map(Fraction, coefficients))
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def add_polynomials(first, second):
    return polynomial([
        (first[index] if index < len(first) else 0)
        + (second[index] if index < len(second) else 0)
        for index in range(max(len(first), len(second)))
    ])


def multiply_polynomials(first, second):
    values = [Fraction(0)] * (len(first) + len(second) - 1)
    for first_index, first_value in enumerate(first):
        for second_index, second_value in enumerate(second):
            values[first_index + second_index] += first_value * second_value
    return polynomial(values)


def sum_polynomials(values):
    total = (Fraction(0),)
    for value in values:
        total = add_polynomials(total, value)
    return total


def interval_bounds(coefficients, upper):
    terms = [value * upper ** index for index, value in enumerate(coefficients) if index > 0]
    return (coefficients[0] + sum(min(0, value) for value in terms),
            coefficients[0] + sum(max(0, value) for value in terms))


def family_certificate():
    half = Fraction(1, 2)
    upper = Fraction(1, 100)
    vectors = {
        'p': ((half, -2), (half, -2), (0, 2), (0, 2)),
        'q': ((half,), (half, -4), (0, 3), (0, 1)),
    }
    expected_coefficients = (
        (0,), (0, 0, -10, 120, -560, 780),
        (0, 0, 5, -360, 1270, -2340),
        (0, 0, -10, 360, -1670, 1560),
        (0, 0, 15, -120, 960),
    )
    masses = {}
    for name, vector in vectors.items():
        require(sum_polynomials(vector) == (1,), 'Family normalization differs')
        for weight in vector:
            lower, _ = interval_bounds(weight, upper)
            require(lower > 0 or (len(weight) == 2 and weight[0] == 0 and weight[1] > 0), 'Family positivity fails')
        for first, second in zip(vector, vector[1:]):
            difference = add_polynomials(first, polynomial([-value for value in second]))
            require(interval_bounds(difference, upper)[0] >= 0, 'Family decreasing order fails')
        masses[name] = [(Fraction(0),)] * 5
        for word in product(range(4), repeat=5):
            weight = (Fraction(1),)
            for label in word:
                weight = multiply_polynomials(weight, vector[label])
            count = len(set(word))
            masses[name][count] = add_polynomials(masses[name][count], weight)
    for length in range(1, 5):
        prefix_difference = add_polynomials(
            sum_polynomials(vectors['q'][:length]),
            polynomial([-value for value in sum_polynomials(vectors['p'][:length])]),
        )
        require(interval_bounds(prefix_difference, upper)[0] >= 0, 'Family majorization fails')
    differences = [add_polynomials(first, polynomial([-value for value in second]))
                   for first, second in zip(masses['p'], masses['q'])]
    require(tuple(differences) == expected_coefficients, 'Family polynomial coefficients differ')
    require(sum_polynomials(differences) == (0,), 'Family total difference is nonzero')
    enclosures = {}
    for count in range(1, 5):
        coefficients = differences[count]
        require(coefficients[:2] == (0, 0), 'Family factor t^2 is absent')
        lower, higher = interval_bounds(coefficients[2:], upper)
        require(higher < 0 if count in (1, 3) else lower > 0, 'Family sign bound fails')
        require(sum(value * upper ** index for index, value in enumerate(coefficients))
                == Fraction(EXPECTED['p'][count] - EXPECTED['q'][count], DENOMINATOR),
                'Family endpoint differs from the integer witness')
        enclosures[str(count)] = [str(lower), str(higher)]
    return {
        'parameter_interval': '0 < t <= 1/100, for every real t',
        'coordinate_coefficients_in_ascending_order': {
            name: [[str(value) for value in coordinate] for coordinate in vector]
            for name, vector in vectors.items()
        },
        'mass_difference_coefficients_in_ascending_order_k0_through_k4': [
            [str(value) for value in coefficients] for coefficients in differences
        ],
        'mass_difference_divided_by_t_squared_enclosures': enclosures,
        'signs_k1_through_k4': ['negative', 'positive', 'negative', 'positive'],
        'coverage': 'Exact polynomial identities and rational coefficientwise interval bounds on the full real interval',
    }


def all_draw_counts_certificate():
    half = Fraction(1, 2)
    vectors = {
        'p': ((half, -2), (half, -2), (0, 2), (0, 2)),
        'q': ((half, 0), (half, -4), (0, 3), (0, 1)),
    }
    terms = {2: [], 3: []}
    for count, coefficients in ((2, (-3, 1, 0)), (3, (3, -2, 1))):
        for size, coefficient in enumerate(coefficients, start=1):
            for name, vector in vectors.items():
                signed_coefficient = coefficient if name == 'p' else -coefficient
                for subset in combinations(range(4), size):
                    base = sum(vector[index][0] for index in subset)
                    slope = sum(vector[index][1] for index in subset)
                    terms[count].append((signed_coefficient, base, slope))
    derivative_bounds = {}
    second_weights = {}
    for count, expression in terms.items():
        derivative_bounds[count] = sum(abs(coefficient * slope ** 3)
                                       for coefficient, base, slope in expression)
        second_weights[count] = {
            str(base): sum(coefficient * slope ** 2
                           for coefficient, term_base, slope in expression if term_base == base)
            for base in sorted({base for coefficient, base, slope in expression})
        }
        for base in {base for coefficient, base, slope in expression}:
            require(sum(coefficient for coefficient, term_base, slope in expression
                        if term_base == base) == 0, 'Constant base weight differs')
            require(sum(coefficient * slope for coefficient, term_base, slope in expression
                        if term_base == base) == 0, 'Linear base weight differs')
    require(derivative_bounds == {2: 684, 3: 1120}, 'Third-derivative bound differs')
    require(second_weights[2] == {'0': 6, '1/2': 4, '1': 0}, 'Second coefficient at two differs')
    require(second_weights[3] == {'0': -6, '1/2': 8, '1': -2}, 'Second coefficient at three differs')
    exact_controls = []
    for draws in range(5, 21):
        upper = Fraction(1, 100 * draws * 2 ** draws)
        margin = draws * (draws - 1) * Fraction(1, 2 ** draws) * upper ** 2
        for count, expression in terms.items():
            coefficients = [sum(coefficient * comb(draws, degree) * base ** (draws - degree)
                                * slope ** degree for coefficient, base, slope in expression)
                            for degree in range(draws + 1)]
            require(coefficients[:2] == [0, 0], 'Constant or linear term differs')
            expected_second = (draws * (draws - 1) * Fraction(1, 2 ** (draws - 3))
                               if count == 2 else draws * (draws - 1)
                               * (Fraction(1, 2 ** (draws - 4)) - 1))
            require(coefficients[2] == expected_second, 'Quadratic coefficient differs')
            actual = sum(coefficient * upper ** degree for degree, coefficient in enumerate(coefficients))
            require(actual >= margin if count == 2 else actual <= -margin, 'Uniform sign margin fails')
        exact_controls.append({'draws': draws, 'parameter': '1/(100*n*2^n)',
                               'common_margin': 'n*(n-1)*2^(-n)*t^2'})
    require(Fraction(114, 800) < 1, 'Positive-sign comparison fails')
    require(Fraction(560, 3 * 100 * 32) == Fraction(7, 120), 'Negative-sign bound differs')
    extension_controls = []
    for draws in (5, 6, 10):
        parameter = Fraction(1, 100 * draws * 2 ** draws)
        margin = draws * (draws - 1) * Fraction(1, 2 ** draws) * parameter ** 2
        for labels in (5, 6, 8):
            extra_labels = labels - 4
            extra_mass = min(margin / (4 * draws), extra_labels * parameter / 4, Fraction(1, 4))
            extended = {}
            for name, vector in vectors.items():
                original = tuple(base + slope * parameter for base, slope in vector)
                extended[name] = tuple((1 - extra_mass) * value for value in original) + (
                    extra_mass / extra_labels,) * extra_labels
                require(sum(extended[name]) == 1 and min(extended[name]) > 0, 'Extension simplex fails')
                require(tuple(sorted(extended[name], reverse=True)) == extended[name], 'Extension order fails')
            require(all(sum(extended['p'][:length]) <= sum(extended['q'][:length])
                        for length in range(1, labels)), 'Extension majorization fails')
            probabilities = {}
            for name, vector in extended.items():
                subset_sums = [sum(sum(vector[index] for index in subset) ** draws
                                   for subset in combinations(range(labels), size))
                               for size in range(1, 4)]
                probabilities[name] = (
                    subset_sums[1] - (labels - 1) * subset_sums[0],
                    subset_sums[2] - (labels - 2) * subset_sums[1]
                    + comb(labels - 1, 2) * subset_sums[0],
                )
            require(probabilities['p'][0] - probabilities['q'][0] >= margin / 2,
                    'Extended positive mass margin fails')
            require(probabilities['p'][1] - probabilities['q'][1] <= -margin / 2,
                    'Extended negative mass margin fails')
            extension_controls.append({'draws': draws, 'labels': labels,
                                       'extra_mass': 'min(margin/(4*n),(N-4)*t/4,1/4)'})
    return {
        'domain': 'Every integer n >= 5 and every real 0 < t <= 1/(100*n*2^n)',
        'quadratic_coefficients': {
            'k2': 'n*(n-1)*2^(3-n)',
            'k3': 'n*(n-1)*(2^(4-n)-1)',
        },
        'second_coefficient_base_weights': second_weights,
        'third_derivative_absolute_slope_cube_sums': derivative_bounds,
        'third_order_remainder_constants': {'k2': '114', 'k3': '560/3'},
        'common_positive_margin': 'n*(n-1)*2^(-n)*t^2',
        'all_dimensions': 'For N > 4 put M=N-4, delta=min(margin/(4*n),M*t/4,1/4); scale the four coordinates by 1-delta and append M coordinates delta/M. Each mass difference changes by at most 2*n*delta <= margin/2.',
        'exact_endpoint_controls': exact_controls,
        'exact_extension_controls': extension_controls,
        'coverage': 'Exact coefficient weights and derivative-bound constants in the written all-n Taylor proof; endpoint controls for n=5,...,20 and nine enlarged-vector pairs are finite checks, not the all-n or all-N proofs.',
    }


def main():
    prefix_sums = {}
    subset_power_sums = {}
    for name, vector in VECTORS.items():
        require(sum(vector) == 100 and min(vector) > 0, 'Invalid probability vector')
        require(tuple(sorted(vector, reverse=True)) == vector, 'Vector order differs')
        prefix_sums[name] = [sum(vector[:length]) for length in range(1, 5)]
        for method in (enumerate_masses, inclusion_exclusion, multinomial_masses):
            require(method(vector) == EXPECTED[name], 'Mass calculation differs: ' + method.__name__)
        require(sum(EXPECTED[name]) == DENOMINATOR, 'Total mass differs from one')
        subset_power_sums[name] = [
            sum(sum(vector[label] for label in subset) ** 5 for subset in combinations(range(4), size))
            for size in range(1, 4)
        ]
    require(all(left <= right for left, right in zip(prefix_sums['p'], prefix_sums['q'])), 'Majorization fails')
    differences = [left - right for left, right in zip(EXPECTED['p'], EXPECTED['q'])]
    require(differences[2] > 0 and differences[3] < 0, 'Threshold obstruction fails')
    cumulative = [sum(differences[:length + 1]) for length in range(5)]
    require(all(value <= 0 for value in cumulative), 'Cumulative order fails')
    result = {
        'schema_version': 1,
        'arithmetic': 'Exact integer witness numerators and rational polynomial arithmetic; no floating-point rounding',
        'draws': 5,
        'labels': 4,
        'denominator': DENOMINATOR,
        'vector_numerators': VECTORS,
        'prefix_sum_numerators': prefix_sums,
        'mass_numerators_k0_through_k4': EXPECTED,
        'subset_power_sum_numerators_j1_through_j3': subset_power_sums,
        'difference_numerators_k0_through_k4': differences,
        'cumulative_difference_numerators_k0_through_k4': cumulative,
        'methods': ['weighted enumeration', 'subset inclusion-exclusion', 'multinomial count enumeration'],
        'words_per_vector': 4 ** 5,
        'occupancy_vectors_per_vector': comb(8, 3),
        'claim': 'The two mass inequalities at counts 2 and 3 exclude every real single-crossing threshold.',
        'coverage': 'All five masses, all decreasing-prefix sums, normalization, positivity and cumulative order for these two vectors.',
        'checks': 'pass',
        'family': family_certificate(),
        'all_draw_counts': all_draw_counts_certificate(),
    }
    output = json.dumps(result, indent=2, sort_keys=True) + '\n'
    destination = Path(__file__).resolve().parent / 'verification.json'
    destination.write_text(output)
    print(output, end='')
    print('verification.json SHA256 ' + hashlib.sha256(output.encode()).hexdigest())


if __name__ == '__main__':
    main()
