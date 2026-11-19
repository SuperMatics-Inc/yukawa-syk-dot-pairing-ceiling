import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

-- A single signed crossing has norm square (x+y)^2. The displayed
-- block consists of the crossing and one copy on each monotone half.
theorem folded_shift_block (x y z u : ℚ) :
    (x - (-y))^2 + (z - u)^2 + ((-z) - (-u))^2 =
      (x+y)^2 + 2*(z-u)^2 := by ring

-- Unit supported h gives crossing x=y=1 and one nonzero same-sign
-- difference. The independently evaluated full-line bubble is 2-(-1).
theorem equal_flavor_transfer_anchor :
    (4 + 2 + (2*3/2)*((1+1)^2 + 2*(1-0)^2) : ℚ) =
      4 + 2 + (2*3)*(2-(-1)) := by norm_num

theorem weak_width_guard :
    (2 : ℚ) * (4046104720/1000000000) * (15/32)^2 < 25/13 := by norm_num

theorem weak_hankel_guard :
    (1 : ℚ) + 1/36 + 1/225 + 16/(49*81) < 26/25 := by norm_num
