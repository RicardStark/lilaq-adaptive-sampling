#import "@local/lilaq:0.6.0" as lq

#set page(
  width: auto,
  height: auto,
  margin: 1pt,
)

#lq.diagram(
  lq.plot-function(
    (-15, 15),
    x => x/calc.sin(x),
    sampling-method: "adaptive",
    mark: auto
  ),
  ylim: (-20, 20)
)