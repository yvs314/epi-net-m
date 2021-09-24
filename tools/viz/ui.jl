"""
Authors: Kara Ignatenko, 2021

Contains useful UI elements for Pluto.jl visualizations.

ui.jl
2021-09-22 v.0.1: First modular version
2021-09-24 v.0.2: Add support for arbitrary steps in the day picker
"""
module VizUI

using HypertextLiteral
using VegaLite

export dayPicker, vegaEmbed

"""
Create an interactive day picker widget for Pluto.jl.
The widget accepts integer values between `MIN` and `MAX`.
In addition to the textbox, buttons are provided to jump to the
first and last values, as well as an arbitrary number of steps
forwards or backwards as given in `steps`.
"""
function dayPicker(min::Int, max::Int, steps::Vector{Int}=[1])
    # create arrays of IDs for backwards and forwards step buttons
    btnIDs = vcat([ "#prev-" * string(step) for step in reverse(steps) ],
                  [ "#next-" * string(step) for step in steps ])


    @htl("""
    <div id="day-picker">
        <button id="first-btn"> first </button>
        $((
            @htl("""<button id="prev-$step"> -$step </button>""") for step in reverse(steps)
        ))
        <input  id="curr-day"/>
        $((
            @htl("""<button id="next-$step"> +$step </button>""") for step in steps
        ))
        <button id="last-btn"> last </button>

        <script>
        const div = currentScript.parentElement;
        const inputField = div.querySelector("#curr-day");

        // min and max day allowed
        const MIN = $(min);
        const MAX = $(max);

        function setValue(newVal) {
            // only change if new number in range and valid
            if (newVal >= MIN && newVal <= MAX && !isNaN(newVal)) {
                div.value = newVal;
                div.dispatchEvent(new CustomEvent("input"));
            }
            // no matter if the value changed, update the textbox
            // to clean up any invalid characters
            inputField.value = div.value;
        }

        setValue(MIN);

        // update value only when Return is pressed or focus is lost
        inputField.addEventListener("change", e => {
            const newVal = parseInt(e.target.value);
            setValue(newVal);
        });

        /*** buttons ***/

        ["#first-btn", "#last-btn"].forEach((s, i) => {
            div.querySelector(s).addEventListener("click", e => {
                const newVal = i == 0 ? MIN : MAX;
                setValue(newVal);
            });
        });

        $(btnIDs).forEach((s, i) => {
            div.querySelector(s).addEventListener("click", e => {
                const newVal = div.value + parseInt(e.target.innerText);
                setValue(newVal);
            });
        });
        </script>
    </div> 
    """)
end

"""
Create the HTML for a Vega-Embed instance that will render
the VegaLite spec `spec`. 
"""
function vegaEmbed(spec::VegaLite.VLSpec)
    # transform spec into JSON
    json = sprint(VegaLite.our_json_print, spec)

    @htl("""
    <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>

    <div id="vis"></div>

    <script>
        const spec = JSON.parse($json);
    
        vegaEmbed("#vis", spec)
            // result.view provides access to the Vega View API
          .then(result => console.log(result))
          .catch(console.warn);
    </script>

    """)

end

end