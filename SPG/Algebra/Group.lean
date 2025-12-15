import SPG.Algebra.Basic

namespace SPG.Algebra

/--
Generate the group closure from a list of generators.
This uses a fixed-point iteration: repeatedly adding all pairwise products
until the set size stabilizes.
-/
partial def generate_group (gens : List SPGElement) : List SPGElement :=
  let rec loop (current : List SPGElement) : List SPGElement :=
    let new_elements := (
      current.flatMap fun g1 =>
      current.flatMap fun g2 =>
      [g1 * g2]
    ).eraseDups
    let combined := (current ++ new_elements).eraseDups
    if combined.length == current.length then
      current
    else
      loop combined
  loop gens

end SPG.Algebra
