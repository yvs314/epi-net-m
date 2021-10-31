#using Revise
# cd("tools/viz")
using Mux, HTTP, DataFrames, FromFile, VegaLite, JSON3
@from "../Oboe/Oboe.jl" using Oboe
@from "./network.jl" using Network

const DIR = @__DIR__

"""Delegate to `app` and respond with the output converted to JSON."""
function jsonify(app, req)
  data = app(req)

  Dict(
    :body => if typeof(data) == VegaLite.VLSpec
      sprint(VegaLite.our_json_print, data)
    else
      JSON3.write(data)
    end,
    :headers => Dict("Content-Type" => "application/json")
  )
end

@app dashboard = (
  Mux.defaults,
  page("/", req -> Mux.fileresponse("$DIR/static/index.html")),
  route("/out/", files("../../out")),
  route("/api/network", mux(jsonify, netHandler)),
  files("$DIR/static", false),
  Mux.notfound()
)
## the default port is 8000
println("ready to serve")
wait(serve(dashboard))