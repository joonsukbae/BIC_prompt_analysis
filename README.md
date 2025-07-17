```mermaid
flowchart LR
  A["<i><code>Event Building</code></i>"] -->|".root tree"| B["<i><code>Event Selection</code></i> <br/> (PID or else?)"]
  B -->|".root tree"| C["<i><code>Calibration Constants</code></i>"]
  B'["DAQ mapping"]-->|".txt"| C
  H["sim results"]-->|".root"| C
  H'["sim mapping"]-->|".txt"| C
  C -->|".root histograms"| D["<i><code>QA Draw Macro</code></i>"]
  C -->|".csv constants"| E["<i><code>Energy Calibration</code></i>"]
  F["new selected data"] -->|".root tree"| E
  E-->G[?]