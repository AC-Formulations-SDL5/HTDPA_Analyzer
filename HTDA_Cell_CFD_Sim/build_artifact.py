"""Builds the self-contained interactive HTML artifact (frame scrubber) from
the cached simulation results, embedding each frame as a base64 PNG."""
import json
import os

import plot_reciprocating_flow as pf

HERE = os.path.dirname(__file__)
OUT_PATH = os.path.join(HERE, "reciprocating_flow_viewer.html")

HTML_TEMPLATE = """<title>Reciprocating Flow — Interactive Viewer</title>
<style>
:root {
  --bg: #ffffff;
  --panel: #f4f5f7;
  --text: #1c1f26;
  --muted: #5b6270;
  --accent: #2f6fed;
  --border: #d9dce2;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #14161b;
    --panel: #1d2028;
    --text: #e7e9ee;
    --muted: #9aa1b0;
    --accent: #5b9dff;
    --border: #30343e;
  }
}
:root[data-theme="dark"] {
  --bg: #14161b;
  --panel: #1d2028;
  --text: #e7e9ee;
  --muted: #9aa1b0;
  --accent: #5b9dff;
  --border: #30343e;
}
* { box-sizing: border-box; }
body {
  background: var(--bg);
  color: var(--text);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  margin: 0;
  padding: 24px 16px 40px;
}
.wrap { max-width: 880px; margin: 0 auto; }
h1 { font-size: 1.25rem; margin: 0 0 4px; }
p.sub { color: var(--muted); margin: 0 0 20px; font-size: 0.92rem; line-height: 1.5; }
.frame-box {
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 12px;
  padding: 12px;
  display: flex;
  justify-content: center;
}
.frame-box img { max-width: 100%; height: auto; border-radius: 6px; display: block; }
.controls {
  display: flex;
  align-items: center;
  gap: 12px;
  margin-top: 16px;
  flex-wrap: wrap;
}
button {
  background: var(--accent);
  color: white;
  border: none;
  border-radius: 8px;
  padding: 8px 16px;
  font-size: 0.9rem;
  cursor: pointer;
}
button:hover { opacity: 0.9; }
input[type="range"] { flex: 1; min-width: 180px; accent-color: var(--accent); }
.time-readout {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  color: var(--text);
  min-width: 92px;
  text-align: right;
  font-size: 0.9rem;
}
.legend {
  margin-top: 14px;
  font-size: 0.82rem;
  color: var(--muted);
  line-height: 1.5;
}
.params {
  margin-top: 18px;
  padding-top: 14px;
  border-top: 1px solid var(--border);
  font-size: 0.82rem;
  color: var(--muted);
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
  gap: 6px 18px;
}
.params b { color: var(--text); }
</style>

<div class="wrap">
  <h1>Reciprocating flow in the rectangular chamber</h1>
  <p class="sub">
    Transient 2D incompressible Navier&ndash;Stokes (NVIDIA Warp <code>warp.fem</code>, Taylor&ndash;Hood Q2&ndash;Q1,
    masked Grid2D + hard velocity Dirichlet BCs). Drag the slider or press Play to scrub through one full
    reciprocation cycle. Color scales are recomputed per frame so instantaneous structure stays visible.
  </p>

  <div class="frame-box">
    <img id="frame-img" src="" alt="flow field frame" />
  </div>

  <div class="controls">
    <button id="play-btn">▶ Play</button>
    <input id="slider" type="range" min="0" max="__NFRAMES_MINUS_1__" value="0" step="1" />
    <span class="time-readout" id="time-readout">t = 0.00 s</span>
  </div>

  <div class="legend">
    Panels: velocity magnitude |u| with direction arrows (mm/s) &middot; kinematic pressure p/&rho; (mm&sup2;/s&sup2;)
    &middot; out-of-plane vorticity (1/s). All fields are cropped to the rectangle interior (pipe stubs excluded).
  </div>

  <div class="params">
    <div><b>Q amplitude</b>: __Q_AMP__ mL/min</div>
    <div><b>Period T</b>: __PERIOD__ s</div>
    <div><b>Rectangle</b>: __RECT_W__ &times; __RECT_H__ mm</div>
    <div><b>Pipe ID</b>: __PIPE_D__ mm</div>
    <div><b>Assumed depth</b>: __DEPTH__ mm</div>
    <div><b>Fluid</b>: water (&nu; = __NU__ mm&sup2;/s)</div>
  </div>
</div>

<script>
const FRAMES = __FRAMES_JSON__;
const TIMES = __TIMES_JSON__;

const img = document.getElementById('frame-img');
const slider = document.getElementById('slider');
const readout = document.getElementById('time-readout');
const playBtn = document.getElementById('play-btn');

let playing = false;
let timer = null;

function showFrame(i) {
  img.src = "data:image/png;base64," + FRAMES[i];
  readout.textContent = "t = " + TIMES[i].toFixed(2) + " s";
}

slider.addEventListener('input', () => showFrame(parseInt(slider.value, 10)));

function stepPlay() {
  let v = parseInt(slider.value, 10) + 1;
  if (v > parseInt(slider.max, 10)) v = 0;
  slider.value = v;
  showFrame(v);
}

playBtn.addEventListener('click', () => {
  playing = !playing;
  if (playing) {
    playBtn.textContent = '⏸ Pause';
    timer = setInterval(stepPlay, 90);
  } else {
    playBtn.textContent = '▶ Play';
    clearInterval(timer);
  }
});

showFrame(0);
</script>
"""


def main():
    d = pf.load_results()
    images, times = pf.build_frame_images(d, max_frames=None)

    html = HTML_TEMPLATE
    html = html.replace("__NFRAMES_MINUS_1__", str(len(images) - 1))
    html = html.replace("__FRAMES_JSON__", json.dumps(images))
    html = html.replace("__TIMES_JSON__", json.dumps([round(x, 4) for x in times]))
    html = html.replace("__Q_AMP__", f"{float(d['q_amp_ml_min']):.0f}")
    html = html.replace("__PERIOD__", f"{float(d['period']):.0f}")
    html = html.replace("__RECT_W__", f"{float(d['rect_w']):.2f}")
    html = html.replace("__RECT_H__", f"{float(d['rect_h']):.1f}")
    html = html.replace("__PIPE_D__", f"{float(d['pipe_d']):.1f}")
    html = html.replace("__DEPTH__", f"{float(d['depth']):.1f}")
    html = html.replace("__NU__", f"{float(d['nu']):.2f}")

    with open(OUT_PATH, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"wrote {OUT_PATH} ({os.path.getsize(OUT_PATH)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
