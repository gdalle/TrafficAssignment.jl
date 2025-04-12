# TrafficAssignment.jl

[![Build Status](https://github.com/gdalle/TrafficAssignment.jl/actions/workflows/Test.yml/badge.svg?branch=master)](https://github.com/gdalle/TrafficAssignment.jl/actions/workflows/Test.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/gdalle/TrafficAssignment.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/gdalle/TrafficAssignment.jl)
[![Dev Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://gdalle.github.io/TrafficAssignment.jl/dev/)

This is a Julia package for studying [traffic assignment](https://en.wikipedia.org/wiki/Route_assignment) on road networks.

## Getting started

To install the latest development version, run this in your Julia Pkg REPL:

```julia
pkg> add https://github.com/gdalle/TrafficAssignment.jl
```

You can then load networks from the [`TransportationNetworks` repository](https://github.com/bstabler/TransportationNetworks):

```jldoctest readme
using TrafficAssignment
ta_data = load_ta_network("SiouxFalls")
ta_data.number_of_zones

# output

24
```

And then you can solve the equilibrium problem:

```jldoctest readme
link_flow, link_travel_time, objective = ta_frank_wolfe(ta_data, log="off", tol=1e-2)
```

The total system travel time can be simply computed as:

```jldoctest readme
using LinearAlgebra
system_travel_time = dot(link_travel_time, link_flow)

# output

7.467895288156025e6
```

## Credits

This package was originally written and maintained by [Changhyun Kwon](http://www.chkwon.net).
