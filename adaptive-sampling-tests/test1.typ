#import "@local/lilaq:0.6.0" as lq

#set page(
  width: auto,
  height: auto,
  margin: 1pt,
)

#lq.diagram(
  lq.plot-function(
    (-15, 15),
    x => calc.sin(x)/x,
    sampling-method: "adaptive",
    mark: auto
  ),
  ylim: (-0.5, 1.2)
)