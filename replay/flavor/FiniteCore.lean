import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

-- Full-line overlap anchor for arbitrary flavor coefficient.
theorem flavor_transfer_anchor (r : ℚ) :
    4 + 2 + (r * 3 / 2) * ((1 + 1)^2 + 2 * (1 - 0)^2) =
      4 + 2 + r * 3 * (2 - (-1)) := by ring

theorem uniform_weak_width_guard :
    (101 / 50 : ℚ) * (4046104720 / 1000000000) * (15 / 32)^2 < 25 / 13 := by
  norm_num

theorem uniform_positive_residual_guard :
    let e : ℚ := 1 / 180
    let c : ℚ := 3 / 4 * (5 / 9 * (16 / 5) - 1 / (101 / 50))
    (c + 15 / 2 * e) / (1 - e * (c + 15 / 2 * e)) < 101 / 100 := by
  norm_num

theorem uniform_first_transfer_guard :
    let e : ℚ := 1 / 180
    (1 + 2 * (99 / 50) * (1 - e / 2)^2) /
        (1 + (1 + 2 * (101 / 50)) * (4 / 5) * e) > 19 / 4 := by
  norm_num
