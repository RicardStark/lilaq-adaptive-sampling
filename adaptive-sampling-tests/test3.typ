#import "@local/lilaq:0.6.0" as lq

#set page(
  width: auto,
  height: auto,
  margin: 1pt,
)

#lq.diagram(
  lq.plot-function(
    (-15, 15),
    x => calc.exp(-2*x)/(x - 1),
    sampling-method: "adaptive",
    mark: auto
  ),
  ylim: (-10, 10)
)