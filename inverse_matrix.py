"""역행렬 계산 프로그램

n x n 정방행렬을 입력받아 두 가지 방법으로 역행렬을 구하고 결과를 비교한다.
  1) 행렬식(여인수/수반행렬)을 이용한 방법
  2) 가우스-조던 소거법을 이용한 방법

추가기능
  - 분수(Fraction) 연산으로 오차 없는 정확한 결과 출력
  - 잘못된 입력에 대한 재입력 처리
  - A x A^-1 = I 검산
"""

from fractions import Fraction


# ---------------------------------------------------------------
# 1. 행렬 입력 기능
# ---------------------------------------------------------------
def input_matrix():
    """정수 n과 n x n 행렬을 행 단위로 입력받아 2차원 리스트로 반환한다."""
    while True:
        try:
            n = int(input("정방행렬의 크기 n을 입력하세요: "))
            if n < 1:
                raise ValueError
            break
        except ValueError:
            print("  [입력 오류] 1 이상의 정수를 입력하세요.")

    print(f"{n} x {n} 행렬을 한 행씩, 공백으로 구분하여 입력하세요.")
    matrix = []
    for i in range(n):
        while True:
            try:
                # Fraction은 정수("3"), 소수("0.5"), 분수("1/3") 입력을 모두 처리한다.
                row = [Fraction(x) for x in input(f"  {i + 1}행: ").split()]
            except (ValueError, ZeroDivisionError):
                print("  [입력 오류] 숫자만 입력하세요.")
                continue
            if len(row) != n:
                print(f"  [입력 오류] 원소를 정확히 {n}개 입력하세요.")
                continue
            matrix.append(row)
            break
    return matrix


# ---------------------------------------------------------------
# 2. 행렬식을 이용한 역행렬 계산 기능
# ---------------------------------------------------------------
def minor(matrix, row, col):
    """row행과 col열을 제거한 소행렬을 반환한다."""
    return [r[:col] + r[col + 1:] for i, r in enumerate(matrix) if i != row]


def determinant(matrix):
    """1행에 대한 여인수 전개로 행렬식을 계산한다(재귀)."""
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    det = 0
    for j in range(n):
        det += (-1) ** j * matrix[0][j] * determinant(minor(matrix, 0, j))
    return det


def inverse_by_determinant(matrix):
    """A^-1 = adj(A) / det(A) 공식으로 역행렬을 계산한다."""
    n = len(matrix)
    det = determinant(matrix)
    if det == 0:
        raise ValueError("행렬식이 0이므로 역행렬이 존재하지 않습니다.")
    if n == 1:
        return [[1 / det]]

    # 여인수 행렬 C[i][j] = (-1)^(i+j) * det(M_ij)
    cofactors = [
        [(-1) ** (i + j) * determinant(minor(matrix, i, j)) for j in range(n)]
        for i in range(n)
    ]
    # 수반행렬은 여인수 행렬의 전치이므로 [j][i]로 접근한다.
    return [[cofactors[j][i] / det for j in range(n)] for i in range(n)]


# ---------------------------------------------------------------
# 3. 가우스-조던 소거법을 이용한 역행렬 계산 기능
# ---------------------------------------------------------------
def inverse_by_gauss_jordan(matrix):
    """첨가행렬 [A | I]를 [I | A^-1]로 변환하여 역행렬을 계산한다."""
    n = len(matrix)
    # 첨가행렬 [A | I] 생성
    aug = [
        list(matrix[i]) + [Fraction(int(i == j)) for j in range(n)]
        for i in range(n)
    ]

    for col in range(n):
        # 피벗이 0이 아닌 행을 찾아 교환한다.
        pivot_row = next((r for r in range(col, n) if aug[r][col] != 0), None)
        if pivot_row is None:
            raise ValueError("피벗을 찾을 수 없으므로 역행렬이 존재하지 않습니다.")
        aug[col], aug[pivot_row] = aug[pivot_row], aug[col]

        # 피벗을 1로 만든다.
        pivot = aug[col][col]
        aug[col] = [x / pivot for x in aug[col]]

        # 피벗 열의 나머지 원소를 0으로 만든다.
        for r in range(n):
            if r != col and aug[r][col] != 0:
                factor = aug[r][col]
                aug[r] = [x - factor * y for x, y in zip(aug[r], aug[col])]

    # 오른쪽 절반이 역행렬이다.
    return [row[n:] for row in aug]


# ---------------------------------------------------------------
# 4. 결과 출력 및 비교 기능
# ---------------------------------------------------------------
def print_matrix(matrix):
    """행렬의 각 열을 정렬하여 출력한다. (Fraction은 3/4 같은 분수 형태로 출력)"""
    cells = [[str(x) for x in row] for row in matrix]
    width = max(len(c) for row in cells for c in row)
    for row in cells:
        print("  [ " + "  ".join(c.rjust(width) for c in row) + " ]")


def multiply(a, b):
    """두 정방행렬의 곱을 반환한다. (검산용)"""
    n = len(a)
    return [
        [sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)]
        for i in range(n)
    ]


def is_identity(matrix):
    """단위행렬 여부를 반환한다."""
    n = len(matrix)
    return all(matrix[i][j] == int(i == j) for i in range(n) for j in range(n))


def main():
    print("=" * 50)
    print(" 역행렬 계산 프로그램")
    print("=" * 50)

    matrix = input_matrix()
    print("\n[입력된 행렬 A]")
    print_matrix(matrix)

    print("\n[방법 1] 행렬식을 이용한 역행렬")
    print(f"  det(A) = {determinant(matrix)}")
    try:
        inv_det = inverse_by_determinant(matrix)
        print_matrix(inv_det)
    except ValueError as e:
        inv_det = None
        print(f"  [오류] {e}")

    print("\n[방법 2] 가우스-조던 소거법을 이용한 역행렬")
    try:
        inv_gj = inverse_by_gauss_jordan(matrix)
        print_matrix(inv_gj)
    except ValueError as e:
        inv_gj = None
        print(f"  [오류] {e}")

    print("\n[결과 비교]")
    if inv_det is None and inv_gj is None:
        print("  두 방법 모두 역행렬이 존재하지 않는다고 판단했습니다. (결과 일치)")
        return
    if inv_det == inv_gj:
        print("  두 방법의 결과가 동일합니다.")
    else:
        print("  두 방법의 결과가 서로 다릅니다.")
        return

    print("\n[검산] A x A^-1")
    product = multiply(matrix, inv_det)
    print_matrix(product)
    if is_identity(product):
        print("  단위행렬이 되므로 올바른 역행렬입니다.")
    else:
        print("  단위행렬이 아닙니다. 계산을 확인하세요.")


if __name__ == "__main__":
    main()
