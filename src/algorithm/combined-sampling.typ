#import "../vec.typ" as vec

// This file implements the combined uniform/adaptive sampling from the paper:
// 
// @article{combinedSampling,
//   author  = "Bayer, Tomáš",
//   title   = "Efficient Plotting the Functions with Discontinuities Based on Combined Sampling",
//   journal = "Geoinformatics FCE CTU",
//   year    = 2018,
//   month   = "Aug",
//   volume  = "17",
//   number  = "2",
//   pages   = "9--30",
//   doi     = "https://doi.org/10.14311/gi.17.2.2",
// }
// 
// For checking against the paper's description, we provide the following dictionary of the paper's symbols and their equivalent in this file:
// 
// Paper        |   Code
// -------------|---------------
// d underbar   |   min-depth
// d overbar    |   max-depth
// s overbar    |   max-splits
// ϵ            |   eps
// α overbar    |   angle-tol
// y overbar    |   inf-tol
// 
// Also, it'll be convenient to define some custom types for explanations, even if Typst doesn't support custom types:
//
// Interval = (float, float)
// 2dPoint = (float, float)
// Polyline = array<2dPoint>
// SingularityResult = (status: "ok", points: Polyline) | (status: "singular", c: float)


/// Helper function for checking for infinite singularities at a point
///
/// -> bool
#let is-extreme-or-nonfinite(
  /// The value to test for being an infinite singularity.
  /// 
  /// -> float
  value,
  /// The lower bound tolerance for the infinite singularity.
  /// 
  /// -> float
  inf-tol
) = {
  return calc.abs(value) > inf-tol or float.is-infinite(value) or float.is-nan(value)
}

/// Unified discontinuity detector. It will check for infinite singularities and call the LR criterion to check for non-smooth behavior.
///
/// -> dictionary
#let detect-singularity(
  /// Function being tested for singularities.
  /// 
  /// -> function: float -> float
  f,
  /// Central point to test for singularity.
  /// 
  /// -> float
  x,
  /// Infinite singularity tolerance.
  /// 
  /// -> float
  inf-tol,
  /// Minimum allowed interval half-width.
  /// 
  /// -> float
  eps,
  /// Lower bound tolerance for LR condition.
  /// 
  /// -> float
  lr-tol,
) = {
  // We make an empty array to store the different stencil values f(x - eps), f(x - eps/2), ... to avoid recomputation. It starts as empty since they will be computed in the infinite singularity test loop, and will only be relevant if no infinite singularity is found
  let f-values = ()

  // We begin checking for infinite singularities
  for k in (-2, -1, 0, 1, 2) {
    let v = f(x + k * eps/2)
    if is-extreme-or-nonfinite(v, inf-tol) {
      return (status: "singular", c: x)
    }
      
    // If no singularity is found, store v to the f-values array for later use
    f-values.push(v)
  }

  ////////////////////////////////
  // Here, we now compute the LR score and test the LR criterion
  ////////////////////////////////
  
  // We compute the local terms
  // f_(r) = 3 f(x) - 4 f(x + eps/2) + f(x + eps) and
  // f_(l) = 3 f(x) - 4 f(x + eps/2) + f(x + eps)
  let f-r = 3*f-values.at(2) - 4*f-values.at(3) + f-values.at(4)
  let f-l = 3*f-values.at(2) - 4*f-values.at(1) + f-values.at(0)

  // and then the numerator |f_r^2 - f_l^2|
  let numer = calc.abs(calc.pow(f-r, 2) - calc.pow(f-l, 2))

  // and denominator |f_r^2 + f_l^2 + 1e-4|
  // Annoyingly, the paper doesn't actually explain it's using the 1e-4 as a stabilizer for 0 denominator cases.
  let denom = calc.abs(calc.pow(f-r, 2) + calc.pow(f-l, 2) + 1e-4)

  // If LR = numer/denom > lr-tol, then the function is probably not smooth at x
  if numer/denom > lr-tol {
    return (status: "singular", c: x)
  } else {
    return (status: "ok", c: x, fc: f-values.at(2))
  }
}

/// Handles a singular interval by locating the discontinuity and either
/// 1. shifting the bounds of the interval if the singularity is close to the edges, or
/// 2. splitting the interval at the singularity if it is in the interior of the interval.
/// 
/// This is Algorithm 4, on page 20, of the original paper.
/// 
/// -> array
#let process-int(
  /// The initial interval to process.
  /// 
  /// -> array
  interval,
  /// The location of the singularity.
  /// 
  /// -> float
  c,
  /// Minimum allowed interval half-width.
  /// 
  /// -> float
  eps,
) = {
  let (a, b) = interval

  // 1. Invalid interval
  if a > b {
    return ()
  }

  // 2. "Empty" interval: if the interval is too small, below tolerance, it is simply dropped
  if calc.abs(b - a) < eps {
    return ()
  }

  // 3. Singular point near the lower bound: shift the lower bound upward
  if a <= c and calc.abs(c - a) <= eps {
    return ((a + eps, b),)
  }

  // 4. Singular point near the upper bound: shift the upper bound downward
  if c <= b and calc.abs(c - b) <= eps {
    return ((a, b - eps),)
  }

  // 5. Singular point strictly inside the interval: split into two children
  if a < c and c < b and calc.abs(c - a) > eps and calc.abs(c - b) > eps {
    let left = (a, c - eps)
    let right = (c + eps, b)

    return (left, right)
  }
}

/// Geometric refinement criterion for accepting a segment.
///
/// Returns true if the polygonal approximation is sufficiently rough compared to the angular tolerance angle-tol to need more refinement.
///
/// -> bool
#let refinement-criterion(
  /// The triplet of points to be tested.
  /// 
  /// -> array
  points,
  /// Maximum allowed angular deviation.
  /// 
  /// -> float
  angle-tol,
) = {
  // Paper states acos ± π, but the latter is only needed if one wants to know rotation direction, which isn't important here, so we just set it to acos
  let p0 = points.at(0)
  let p1 = points.at(1)
  let p2 = points.at(2)

  let u = vec.subtract(p2, p1)
  let v = vec.subtract(p0, p1)

  let dot = vec.inner(u, v)
  let norm-u = calc.sqrt(vec.inner(u, u))
  let norm-v = calc.sqrt(vec.inner(v, v))

  // Protect against 0 denominator
  if norm-u * norm-v <= 0 {return false}

  let cosine = calc.abs(dot) / (norm-u * norm-v)
  // Protect against floating point values going beyond [-1, 1]
  let cosine = if cosine > 1 { 1 } else if cosine < -1 { -1 } else { cosine }

  return calc.acos(cosine) > angle-tol
}

/// Recursive combined-sampling routine used inside cs-init. This is the paper's Algorithm 2 on page 14.
///
/// -> dictionary
#let cs(
  /// Function being interpolated.
  /// 
  /// -> function
  f,
  /// Interval [a, b] to interpolate the function `f` on.
  /// 
  /// -> array
  interval,
  /// Value of f(a), to avoid recomputation.
  /// 
  /// -> float
  fa,
  /// Value of f(b), to avoid recomputation.
  /// 
  /// -> float
  fb,
  /// Polyline being built.
  /// 
  /// -> array
  polyline,
  /// Current recursion depth.
  /// 
  /// -> int
  depth,
  /// Minimum required recursion depth.
  /// 
  // -> int
  depth-min,
  /// Maximum allowed recursion depth.
  /// 
  /// -> int
  depth-max,
  /// Minimum allowed interval half-width.
  /// 
  // -> float
  eps,
  /// Lower bound tolerance for LR condition.
  /// 
  /// -> float
  lr-tol,
  /// Infinite singularity tolerance.
  /// 
  /// -> float
  inf-tol,
  /// Maximum allowed angular deviation.
  /// 
  /// -> float
  angle-tol,
  /// Seed to use for random number generation. If set to auto, the seed will be computed as the hash of the current state of the polyline `polyline`.
  /// 
  /// -> auto | int
  seed: auto,
) = {
  // Use suiji for random number generation
  import "@preview/suiji:0.5.1"

  if seed == auto {
    seed = polyline
      .enumerate()
      .map(((index, (l, r))) => array(((l + r)).to-bytes()).sum() * index)
      .sum()
  }

  let rng = suiji.gen-rng-f(seed)

  let (a, b) = interval

  // If either the recursion depth is exceeded or the current interval is too small, end the process early
  if depth > depth-max or b - a < eps {
    return (status: "ok", points: polyline)
  }

  let (r0, r1, r2) = suiji.uniform-f(rng, low: 0.45, high: 0.55, size: 3).at(1)

  let x0 = a + 0.5*r0*(b - a)

  let x0-state = detect-singularity(f, x0, inf-tol, eps, lr-tol)

  if x0-state.status == "singular" {
    return x0-state
  }

  let x1 = a + r1*(b - a)

  let x1-state = detect-singularity(f, x1, inf-tol, eps, lr-tol)

  if x1-state.status == "singular" {
    return x1-state
  }

  let x2 = a + 1.5*r2*(b - a)

  let x2-state = detect-singularity(f, x2, inf-tol, eps, lr-tol)

  if x2-state.status == "singular" {
    return x2-state
  }

  // If no singularities are found, pull y0 = f(x0), y1 = f(x1), and y2 = f(x2) from their status dictionaries, and make the points pa, p0, p1, p2, and pb
  let (pa, p0, p1, p2, pb) = ((a, fa), (x0, x0-state.fc), (x1, x1-state.fc), (x2, x2-state.fc), (b, fb))

  // Compute angular deviation of (pa, p0, p1)
  let alpha0 = refinement-criterion((pa, p0, p1), angle-tol)

  if alpha0 or (depth < depth-min) {
    let child-state = cs(
      f,
      (a, x0),
      fa,
      p0.at(1),
      polyline,
      depth + 1,
      depth-min,
      depth-max,
      eps,
      lr-tol,
      inf-tol,
      angle-tol,
    )

    if child-state.status == "singular" {
      return child-state
    }

    // Update polyline with child-state
    polyline = child-state.points
  }

  polyline.push(p0)

  /// Compute angular deviation of (p0, p1, p2)
  let alpha1 = refinement-criterion((p0, p1, p2), angle-tol)

  if alpha0 or alpha1 or (depth < depth-min) {
    let child-state = cs(
      f,
      (x0, x1),
      p0.at(1),
      p1.at(1),
      polyline,
      depth + 1,
      depth-min,
      depth-max,
      eps,
      lr-tol,
      inf-tol,
      angle-tol,
    )

    if child-state.status == "singular" {
      return child-state
    }

    // Update polyline with child-state
    polyline = child-state.points
  }

  polyline.push(p1)

  /// Compute angular deviation of (p1, p2, pb)
  let alpha2 = refinement-criterion((p1, p2, pb), angle-tol)

  if alpha1 or alpha2 or (depth < depth-min) {
    let child-state = cs(
      f,
      (x1, x2),
      p1.at(1),
      p2.at(1),
      polyline,
      depth + 1,
      depth-min,
      depth-max,
      eps,
      lr-tol,
      inf-tol,
      angle-tol,
    )

    if child-state.status == "singular" {
      return child-state
    }

    // Update polyline with child-state
    polyline = child-state.points
  }

  polyline.push(p2)

  if alpha2 or (depth < depth-min) {
    let child-state = cs(
      f,
      (x2, b),
      p2.at(1),
      fb,
      polyline,
      depth + 1,
      depth-min,
      depth-max,
      eps,
      lr-tol,
      inf-tol,
      angle-tol,
    )

    if child-state.status == "singular" {
      return child-state
    }

    // Update polyline with child-state
    polyline = child-state.points
  }
  
  return (status: "ok", points: polyline)
}

/// Tries to build a polygonal approximation for a function over a single interval.
///
/// This is the per-interval entry point corresponding to the paper's Algorithm 1, on page 13. It samples the interval, checks the geometric criterion, and either accepts the interval or signals a singularity that must be handled.
/// 
/// -> dictionary
#let cs-init(
  /// Function being tested for singularities.
  /// 
  /// -> function
  f,
  /// Interval to interpolate the function `f` on.
  /// 
  /// -> array
  interval,
  /// Minimum required recursion depth.
  /// 
  /// -> int
  depth-min,
  /// Maximum allowed recursion depth.
  /// 
  /// -> int
  depth-max,
  /// Minimum allowed interval half-width.
  /// 
  /// -> float
  eps,
  /// Lower bound tolerance for LR condition.
  /// 
  /// -> float
  lr-tol,
  /// Infinite singularity tolerance.
  /// 
  /// -> float
  inf-tol,
  /// Maximum allowed angular deviation.
  /// 
  /// -> float
  angle-tol,
) = {
  let (a, b) = interval

  // Test if singularity is at or near interval bounds
  let a-state = detect-singularity(f, a, inf-tol, eps, lr-tol)

  if a-state.status == "singular" {
    return a-state
  }

  let b-state = detect-singularity(f, b, inf-tol, eps, lr-tol)

  if b-state.status == "singular" {
    return b-state
  }

  // initialize the Polyline object
  let polyline = ((a, a-state.fc),)

  // Process the interval through combined sampling to either return the full polyline, or another singularity point for processing
  let polyline-state = cs(
    f,
    interval,
    a-state.fc,
    b-state.fc,
    polyline,
    1,
    depth-min,
    depth-max,
    eps,
    lr-tol,
    inf-tol,
    angle-tol,
  )

  // if the combined sampling returns a singularity point, return the singularity state for feeding into process-int
  if polyline-state.status == "singular" {
    return polyline-state
  }

  // if the combined sampling returns a complete polyline, add the endpoint of the interval last
  let polyline-final = polyline-state.points
  polyline-final.push((b, b-state.fc))

  return (status: "ok", points: polyline-final)
}

/// Computes a polygonal approximation for a function over an interval, including infinite singularities and non-smooth behavior. This is the paper's algorithm 3 on page 18. Returns an array of Polylines.
///
/// -> array
#let cs-stack(
  /// Function being interpolated.
  /// 
  /// -> function
  f,
  /// Interval to interpolate the function `f` on.
  /// 
  /// -> array
  interval,
  /// Maximum split depth allowed for handling singularity behavior.
  /// 
  /// -> int
  max-splits: 20,
  /// Minimum required recursion depth.
  /// 
  /// -> int
  depth-min: 0,
  /// Maximum allowed recursion depth.
  /// 
  /// -> int
  depth-max: 6,
  /// Minimum allowed interval half-width.
  /// 
  /// -> float
  eps: 0.0001,
  /// Lower bound tolerance for LR condition.
  /// 
  /// -> float
  lr-tol: 0.8,
  /// Infinite singularity tolerance.
  /// 
  /// -> float
  inf-tol: 1e5,
  /// Maximum allowed angular deviation.
  /// 
  /// -> float
  angle-tol: 1deg,
) = {
  // stack for tracking intervals that must be processed
  let stack = (interval,)

  // final object storing an array of Polyline objects for plotting
  let polyline-list = ()

  // counter for the number of splits, if it goes above the threshold max-splits, dump the entire process as a failure
  let splits = 0

  while stack.len() > 0 {
    if splits > max-splits {
      return ()
    }

    let temp-interval = stack.pop()
    // Try and make a polygonal approximation of f on the popped interval
    let res = cs-init(f, temp-interval, depth-min, depth-max, eps, lr-tol, inf-tol, angle-tol)

    if res.status == "ok" {
      // If the polygonal approximation had no issues, add it to the final colleciton
      polyline-list.push(res.points)
    } else if res.status == "singular" {
      // If a singularity is found in the interval, process the interval by either shifting the interval bounds or splitting it
      let children = process-int(temp-interval, res.c, eps)
      // If the interval is split, iterate the split counter
      if children.len() == 2 {
        splits += 1
      }
      // Push the new interval(s) back onto the stack for processing
      for child in children {
        stack.push(child)
      }
    }
  }

  return polyline-list
}