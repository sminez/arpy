"""
arpy (Absolute Relativity in Python)
Copyright (C) 2016-2018 Innes D. Anderson-Morrison All rights reserved.

Tools for visualising results from calculations and investigations.
"""
import os
import tempfile
import time
import webbrowser

from ..algebra.config import config
from ..algebra.data_types import Alpha, Term
from ..algebra.operations import full
from .lexparse import ARContext

# HTML/JS for js_cayley visualisation
HTML = """\
<head>
  <meta charset="utf-8">
  <script src="https://d3js.org/d3.v4.min.js"></script>
  <script src="https://d3js.org/d3-selection-multi.v0.4.min.js"></script>
  <link href="https://fonts.googleapis.com/css?family=Open+Sans+Condensed:300" rel="stylesheet">

    <style>
        .my-text {
            font-family: 'Open Sans Condensed', sans-serif;
        }

        .cell {
          stroke: #404040;
          shape-rendering: crispEdges;
        }

        div.tooltip {
          position: absolute;
          text-align: center;
          width: 60px;
          height: 35px;
          padding: 2px;
          font: 14px sans-serif;
          background: #EBDBB2;
          border: 1px;
          border-radius: 8px;
          border-style: solid;
          border-color: #404040;
          pointer-events: none;
        }
    </style>

</head>

<body>
    <br>
    <span class="my-text"><H1><u>AR Cayley Table</u></H1>
    <H2>Hover for a summary: click to toggle between value and sign information.</H2></span>
    Algebra components: REPLACE_ALLOWED

    <div id='chart'></div>

    <script>
    d3.functor = function functor(v) {
      return typeof v === "function" ? v : function() {
        return v;
      };
    };

    function getColor(d) {
        if (activeScale == colorScale) {
            return activeScale(d.value);
        } else {
            return activeScale(d.sign);
        }
    }

    var margin = {top: 40, right: 40, bottom: 40, left: 40},
        width = 800 - margin.left - margin.right,
        height = 800 - margin.top - margin.bottom;

    var signColorScale = d3.scaleOrdinal()
      .domain(["+ve", "-ve"])
      .range(["#fbf1c7", "#504945"]);

    var colorScale = d3.scaleOrdinal()
      .domain(["p", "23", "31", "12",
               "0", "023", "031", "012",
               "123", "1", "2", "3",
               "0123", "01", "02", "03"])
      .range(["#dbc78f", "#4a778f", "#6fa0b0", "#8fb5bf",
              "#872e2e", "#76718f", "#9c8ead", "#b49bbf",
              "#4c4861", "#bd4442", "#d66847", "#d98e66",
              "#dbaa56", "#498c6b", "#84b394", "#9fc4a0"]);

    var activeScale = colorScale;

    var combined = REPLACE_COMPONENTS

    var xScale = d3.scaleLinear()
        .range([0, width])
        .domain([0,combined[0].length]);

    var yScale = d3.scaleLinear()
        .range([0, height])
        .domain([0,combined.length]);

    var div = d3.select("body").append("div")
        .attr("class", "tooltip")
        .style("opacity", 0);

    var svg = d3.select("#chart").append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
        .on("click", function(){
            activeScale = (activeScale == colorScale) ? signColorScale : colorScale;
            row.selectAll(".cell").style("fill", function(d) { return getColor(d); });
        });

    var row = svg.selectAll(".row")
        .data(combined)
        .enter().append("svg:g")
        .attr("class", "row");

    var col = row.selectAll(".cell")
        .data(function (d,i) { return d.map(function(a) { return {value: a.val, sign: a.sign, row: i}; } ) })
        .enter().append("svg:rect")
        .attr("class", "cell")
        .attr("x", function(d, i) { return xScale(i); })
        .attr("y", function(d, i) { return yScale(d.row); })
        .attr("width", xScale(1))
        .attr("height", yScale(1))
        .style("fill", function(d) { return getColor(d); })
        .on("mouseover", function(d, i) {
            div.transition()
                .duration(100)
                .style("opacity", .9);
            div.html("α" + d.value + "<br/>" + d.sign)
                .style("left", (d3.event.pageX) + "px")
                .style("top", (d3.event.pageY - 28) + "px");
        })
        .on("mouseout", function(d, i) {
            div.transition()
                .duration(300)
                .style("opacity", 0);
            d3.select(this).style("fill", function() {
                return "" + getColor(d) + "";
            });
        });

    </script>
</body>
"""


def cayley(op=full, padding=6, cfg=config):
    """
    Print current Cayley table to the terminal allowing for specification
    of the operation used to compute the table.

    Any function that accepts two Alphas can be passed as op.
    """
    comps = (
        " ".join(
            [
                str(op(Alpha(a, cfg=cfg), Alpha(b, cfg=cfg), cfg=cfg)).rjust(padding)
                for b in cfg.allowed
            ]
        )
        for a in cfg.allowed
    )
    for comp in comps:
        print(comp)


def sign_cayley(op=full, cfg=config):
    """
    Print +-1 signs for the current Cayley table to the terminal allowing
    for specification of the operation used to compute the table.

    Any function that accepts two Alphas can be passed as op.
    """
    divider = "      " + "".join("+---------" for _ in range(4)) + "+"
    comps = (
        " ".join(
            [
                "■" if op(Alpha(x, cfg=cfg), Alpha(y, cfg=cfg), cfg=cfg).sign == -1 else "□"
                for y in cfg.allowed
            ]
        )
        for x in cfg.allowed
    )

    print("          ", "         ".join(["B", "T", "A", "E"]))
    print(divider)

    for i, comp in enumerate(comps):
        comp = "| ".join(comp[n : n + 8] for n in range(0, len(comp), 8))
        print(str(Alpha(cfg.allowed[i], cfg=cfg)).ljust(5), "|", comp, "|")
        # Divide after each zet
        if (i + 1) % 4 == 0:
            print(divider)


def _4block(rows, cols, op, cfg):
    """Visualise a 4x4 block of 4 elements acting on 4 others"""
    block = []
    for r in rows:
        comps = [op(Alpha(r, cfg=cfg), Alpha(c, cfg=cfg), cfg=cfg).sign for c in cols]
        block_row = " ".join(["□" if c == 1 else "■" for c in comps])
        block.append("|" + block_row + "|")
    return block


def sign_distribution(op=full, cfg=config):
    """
    By calculating one term for each of the 5 grouped components of the
    cayley table, look at how each metric / allowed set of indices affects
    the overall structure of the algebra.
    """
    allowed = cfg.allowed
    bs, xs, zs = allowed[0:16:4], allowed[1:16:4], allowed[3:16:4]

    blocks = []

    row_cols = {"∂e": (bs, bs), "∂Ξ": (bs, xs), "∇": (xs, bs), "∇•": (xs, xs), "∇x": (zs, xs)}

    for name in ["∂e", "∂Ξ", "∇", "∇•", "∇x"]:
        rows, cols = row_cols[name]
        blocks.append((name, _4block(rows, cols, op, cfg)))

    for i in range(4):
        for block in blocks:
            name = block[0].rjust(3) if i == 0 else "   "
            print(name, block[1][i], end=" ")
        print("")


def js_cayley(op=full, cfg=config):
    """
    Compute a full Cayley table for a given config and then
    render the result using an html/js template.

    NOTE: The use of this function requires an internet connection
          in order to load the JavaScript library D3.
    """
    tmp = '{"val": "%s", "sign": "%s"}'

    comps = (
        [op(Alpha(x, cfg=cfg), Alpha(y, cfg=cfg), cfg=cfg) for y in cfg.allowed]
        for x in cfg.allowed
    )

    # This strange format is the JSON structure required by the JS script
    # to parse the data points and generate the Cayley table
    json_comps = [
        "[" + ",".join([tmp % (c._index, "+ve" if c.sign == 1 else "-ve") for c in row]) + "]"
        for row in comps
    ]
    json_comps = "[" + ",\n".join(json_comps) + "];"

    # A horrible hack that really should be using a templating engine but
    # bringing in a dependency just for this seems like overkill...
    html = HTML.replace("REPLACE_ALLOWED", ", ".join(["α{}".format(a) for a in cfg.allowed]))
    html = html.replace("REPLACE_COMPONENTS", json_comps)

    # Write out to a temp file so that we can open it with the browser
    with tempfile.NamedTemporaryFile(mode="w", suffix=".html", delete=False) as f:
        f.write(html)

    webbrowser.open("file://" + f.name)
    time.sleep(1)  # Make sure that the file is opened before removing it
    os.remove(f.name)  # remove the temp file
    print(
        (
            "An internet connection is required to generate the output "
            "due to the use of the D3 JavaScript library."
        )
    )


def op_block(rows, cols, op=full, cfg=config):
    """
    Visualise the sign of component interactions under a binary operation on
    Alphas. In addition to setting `op` to one of the built in binary
    operations, you can also pass your own function with the following
    signature:
        op(Alpha, Alpha, ARConfig) -> Alpha

    Alternatively, you may also pass a valid `ar` string defining how to combine
    each pair of components in terms of i and j:
        "i ^ (j!) ^ a0123"
    """

    def _alpha(x):
        if isinstance(x, Alpha):
            return x
        elif isinstance(x, Term):
            return x.alpha
        else:
            raise ValueError("Must pass a MultiVector or a list of Alphas/Terms")

    def _block_func(s):
        def func(i, j, cfg=config):
            with ARContext(cfg=cfg) as ar:
                return ar(s)

        return func

    if isinstance(op, str):
        # Convert the string to an ar function
        op = _block_func(op)

    block = ""
    for r in rows:
        comps = [op(r.alpha, c.alpha, cfg).sign for c in cols]
        block_row = " ".join(["□" if c == 1 else "■" for c in comps])
        block += "|" + block_row + "|\n"
    print(block)
