# Functions related to the counting bound.


module CountingBound

export countingBound

"""
Counting bound, based on Shannon's argument.
  m: number of edges
  w: number of 'wires' -- that is, log2(number of functions),
    which is the number of bits needed to specify a function
  Returns: average number of NAND gates (with unbounded fan-in)
    required to compute any of those functions.
    (This may not be an integer).
  ??? check this? (e.g. countingBound(3,9) should be 1, I think...
"""
function countingBound(m, w)
  b = m - 0.5
  # the "-1" here is because this is the average, not the max.
  sqrt(2*w + b*b) - b - 1
end




end

