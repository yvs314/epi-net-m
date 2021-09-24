module VizUI

using HypertextLiteral
using VegaLite

export dayPicker, vegaEmbed


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