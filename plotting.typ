#import "@local/lilaq:0.6.0" as lq

#lq.diagram(
  lq.plot-function(
    (-calc.pi, calc.pi),
    x => calc.sin(x),
    sampling-method: "adaptive",
  )
)