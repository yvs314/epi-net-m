module VizUI

using HypertextLiteral
using VegaLite

export dayPicker, vegaEmbed


function dayPicker(min::Int, max::Int)
    @htl("""
    <div id="day-picker">
        <button id="first-btn"> first </button>
        <button id="prev-btn"> prev </button>
        <input  id="curr-day"/>
        <button id="next-btn"> next </button>
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
            }
            // no matter if the value changed, update the Pluto value
            // and the textbox to ensure everything matches
            inputField.value = div.value;
            div.dispatchEvent(new CustomEvent("input"));
        }

        setValue(MIN);

        // update value only when Return is pressed or focus is lost
        inputField.addEventListener("change", e => {
            const newVal = parseInt(e.target.value);
            setValue(newVal);
        });

        /*** buttons ***/

        ["#first-btn", "#prev-btn", "#next-btn", "#last-btn"].forEach((s, i) => {
            div.querySelector(s).addEventListener("click", e => {
                const newVal = i == 0 ? MIN :
                            i == 1 ? div.value - 1 :
                            i == 2 ? div.value + 1 :
                                     MAX;

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