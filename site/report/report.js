// Precipitation intercomparison report — interactive tables and charts.
// Plain JS + SVG. All data comes from ./data/*.json written by scripts/analyze.py.

const STYLE = {
  gpcc_full:       { color: "#1f2324", dash: "" },
  gpcc_monitoring: { color: "#1f2324", dash: "4 3" },
  cru:             { color: "#9aa3a4", dash: "" },
  cpc:             { color: "#9aa3a4", dash: "2 2" },
  precl:           { color: "#9aa3a4", dash: "6 3" },
  udel:            { color: "#9aa3a4", dash: "1 3" },
  chirps_v2:       { color: "#eb6834", dash: "" },
  chirp_v2:        { color: "#eb6834", dash: "5 3" },
  chirps_v3:       { color: "#1baf7a", dash: "" },
  chirp_v3:        { color: "#1baf7a", dash: "5 3" },
  imerg_final:     { color: "#4a3aa7", dash: "" },
  imerg_late:      { color: "#4a3aa7", dash: "5 3" },
  era5:            { color: "#2a78d6", dash: "" },
  asap:            { color: "#e87ba4", dash: "" },
};
const DIV = ["#8f1d1d", "#d03b3b", "#ef9a8a", "#f0efec", "#86b6ef", "#2a78d6", "#104281"];
const SEQ = ["#eef4fc", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"];
const DARK = new Set(["#8f1d1d", "#d03b3b", "#2a78d6", "#104281", "#3987e5", "#256abf", "#184f95", "#0d366b"]);

let META, LABEL = {}, SHORT = {};
const $ = (s, el = document) => el.querySelector(s);
const el = (tag, attrs = {}, ...kids) => {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") e.className = v; else if (k === "html") e.innerHTML = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v); else e.setAttribute(k, v);
  }
  for (const k of kids) if (k != null) e.append(k);
  return e;
};
const svgEl = (tag, attrs = {}) => {
  const e = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  return e;
};
const getJSON = (name) => fetch(`data/${name}.json`).then((r) => { if (!r.ok) throw new Error(name); return r.json(); });
const fmt = (v, d = 2) => (v == null || !isFinite(v) ? "–" : (+v).toFixed(d));
const pct = (v, d = 0) => (v == null || !isFinite(v) ? "–" : `${v > 0 ? "+" : ""}${(+v).toFixed(d)}%`);

// ----------------------------------------------------------------- colour scales
function divBins(v, edges) { // edges ascending, length 6 -> 7 classes
  if (v == null || !isFinite(v)) return null;
  let i = 0; while (i < edges.length && v > edges[i]) i++;
  return DIV[i];
}
const ratioColor = (r) => divBins(r, [0.67, 0.8, 0.95, 1.05, 1.25, 1.5]);
const trendColor = (t) => divBins(t, [-20, -10, -2, 2, 10, 20]);
function seqColor(v, lo = 0, hi = 1) {
  if (v == null || !isFinite(v)) return null;
  const f = Math.max(0, Math.min(0.999, (v - lo) / (hi - lo)));
  return SEQ[Math.floor(f * SEQ.length)];
}
function cellStyle(td, color) {
  if (!color) { td.classList.add("na"); return; }
  td.style.background = color;
  if (DARK.has(color)) td.classList.add("ink-light");
}
function scaleLegend(colors, labels) {
  const s = el("div", { class: "scale" });
  s.append(el("span", {}, labels[0]));
  const bar = el("span", { class: "bar" });
  colors.forEach((c) => bar.append(el("span", { style: `background:${c}` })));
  s.append(bar, el("span", {}, labels[1]));
  return s;
}

// ----------------------------------------------------------------- tables
function heatTable({ rows, cols, value, color, format, rowLabel = (r) => LABEL[r] || r, colLabel = (c) => c, rotate = false, title }) {
  const t = el("table", { class: "t" });
  const hr = el("tr", {}, el("th", {}, title || ""));
  cols.forEach((c) => hr.append(el("th", rotate ? { class: "rot" } : {}, rotate ? el("div", {}, colLabel(c)) : colLabel(c))));
  t.append(el("thead", {}, hr));
  const tb = el("tbody");
  rows.forEach((r) => {
    const tr = el("tr", {}, el("td", { class: "rowh" }, rowLabel(r)));
    cols.forEach((c) => {
      const v = value(r, c);
      const td = el("td", {}, format(v, r, c));
      cellStyle(td, color(v, r, c));
      td.title = `${rowLabel(r)} · ${colLabel(c)}: ${format(v, r, c)}`;
      tr.append(td);
    });
    tb.append(tr);
  });
  t.append(tb);
  return el("div", { class: "tablewrap" }, t);
}
function seg(options, onChange, initial) {
  const s = el("span", { class: "seg", role: "group" });
  options.forEach(([val, label]) => {
    const b = el("button", { type: "button", "aria-pressed": String(val === initial) }, label);
    b.addEventListener("click", () => {
      s.querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", "false"));
      b.setAttribute("aria-pressed", "true");
      onChange(val);
    });
    s.append(b);
  });
  return s;
}
function selectBox(options, onChange, initial) {
  const s = el("select");
  options.forEach(([val, label]) => s.append(el("option", { value: val }, label)));
  if (initial != null) s.value = initial;
  s.addEventListener("change", () => onChange(s.value));
  return s;
}

// ----------------------------------------------------------------- line chart
function lineChart(host, { x, series, yLabel, yFormat = (v) => fmt(v, 0), refY = null, height = 300 }) {
  host.innerHTML = "";
  const W = 820, H = height, m = { l: 54, r: 14, t: 10, b: 28 };
  const shown = new Set(series.map((s) => s.key));
  const legend = el("div", { class: "legend" });
  series.forEach((s) => {
    const b = el("button", { type: "button", "aria-pressed": "true" },
      el("span", { class: "sw", style: `border-top-color:${s.color};border-top-style:${s.dash ? "dashed" : "solid"}` }), s.label);
    b.addEventListener("click", () => {
      if (shown.has(s.key)) shown.delete(s.key); else shown.add(s.key);
      b.setAttribute("aria-pressed", String(shown.has(s.key)));
      draw();
    });
    legend.append(b);
  });
  const box = el("div", { class: "chart" });
  const tip = el("div", { class: "tip" });
  host.append(legend, box);
  box.append(tip);
  function draw() {
    box.querySelectorAll("svg").forEach((n) => n.remove());
    const act = series.filter((s) => shown.has(s.key));
    const vals = act.flatMap((s) => s.values.filter((v) => v != null && isFinite(v)));
    if (refY != null) vals.push(refY);
    if (!vals.length) return;
    let lo = Math.min(...vals), hi = Math.max(...vals);
    const pad = (hi - lo) * 0.08 || Math.abs(hi) * 0.1 || 1; lo -= pad; hi += pad;
    const X = (i) => m.l + (i / (x.length - 1)) * (W - m.l - m.r);
    const Y = (v) => m.t + (1 - (v - lo) / (hi - lo)) * (H - m.t - m.b);
    const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": yLabel });
    const g = svgEl("g", { class: "grid" }), ax = svgEl("g", { class: "axis" });
    const ticks = niceTicks(lo, hi, 5);
    ticks.forEach((tv) => {
      g.append(svgEl("line", { x1: m.l, x2: W - m.r, y1: Y(tv), y2: Y(tv) }));
      const t = svgEl("text", { x: m.l - 6, y: Y(tv) + 4, "text-anchor": "end" }); t.textContent = yFormat(tv); ax.append(t);
    });
    const step = Math.ceil(x.length / 10);
    x.forEach((xv, i) => {
      if (i % step) return;
      const t = svgEl("text", { x: X(i), y: H - 8, "text-anchor": "middle" }); t.textContent = xv; ax.append(t);
    });
    const yl = svgEl("text", { x: 12, y: m.t + (H - m.t - m.b) / 2, transform: `rotate(-90 12 ${m.t + (H - m.t - m.b) / 2})`, "text-anchor": "middle" });
    yl.textContent = yLabel; ax.append(yl);
    svg.append(g, ax);
    if (refY != null) svg.append(svgEl("line", { x1: m.l, x2: W - m.r, y1: Y(refY), y2: Y(refY), stroke: "#5e6a6b", "stroke-width": 1, "stroke-dasharray": "3 3" }));
    act.forEach((s) => {
      let d = "", pen = false;
      s.values.forEach((v, i) => {
        if (v == null || !isFinite(v)) { pen = false; return; }
        d += `${pen ? "L" : "M"}${X(i).toFixed(1)},${Y(v).toFixed(1)}`; pen = true;
      });
      svg.append(svgEl("path", { d, fill: "none", stroke: s.color, "stroke-width": s.key === "gpcc_full" ? 2.4 : 2, "stroke-dasharray": s.dash || "", "stroke-linejoin": "round" }));
    });
    const cross = svgEl("line", { y1: m.t, y2: H - m.b, stroke: "#5e6a6b", "stroke-width": 1, visibility: "hidden" });
    svg.append(cross);
    const hit = svgEl("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent" });
    svg.append(hit);
    hit.addEventListener("mousemove", (ev) => {
      const r = svg.getBoundingClientRect();
      const px = ((ev.clientX - r.left) / r.width) * W;
      const i = Math.max(0, Math.min(x.length - 1, Math.round(((px - m.l) / (W - m.l - m.r)) * (x.length - 1))));
      cross.setAttribute("x1", X(i)); cross.setAttribute("x2", X(i)); cross.setAttribute("visibility", "visible");
      const rows = act.map((s) => ({ s, v: s.values[i] })).filter((o) => o.v != null && isFinite(o.v)).sort((a, b) => b.v - a.v);
      tip.innerHTML = `<div class="tt">${x[i]}</div>` + rows.map((o) =>
        `<div class="row"><span class="k"><i style="background:${o.s.color}"></i>${o.s.label}</span><span>${yFormat(o.v)}</span></div>`).join("");
      tip.style.display = "block";
      const left = (X(i) / W) * r.width;
      tip.style.left = `${left + 14 + 220 > r.width ? left - 230 : left + 14}px`;
      tip.style.top = "8px";
    });
    hit.addEventListener("mouseleave", () => { tip.style.display = "none"; cross.setAttribute("visibility", "hidden"); });
    box.prepend(svg);
  }
  draw();
}
function niceTicks(lo, hi, n) {
  const span = hi - lo, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 2.5, 5, 10].map((k) => k * mag).find((s) => span / s <= n) || 10 * mag;
  const out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}
const seriesFor = (keys, valuesOf) => keys.map((k) => ({ key: k, label: SHORT[k] || k, color: STYLE[k]?.color || "#999", dash: STYLE[k]?.dash || "", values: valuesOf(k) }));

// ----------------------------------------------------------------- sections
function renderProducts(meta) {
  const t = el("table", { class: "t" });
  t.append(el("thead", {}, el("tr", {}, ...["Product", "Type", "Months in this analysis"].map((h) => el("th", {}, h)))));
  const tb = el("tbody");
  meta.products.forEach((p) => tb.append(el("tr", {},
    el("td", { class: "rowh" }, el("span", { class: "sw", style: `display:inline-block;width:16px;border-top:2px ${STYLE[p.key]?.dash ? "dashed" : "solid"} ${STYLE[p.key]?.color};vertical-align:middle;margin-right:6px` }), p.label),
    el("td", { style: "text-align:left" }, p.family), el("td", { style: "text-align:left" }, `${p.span[0]} – ${p.span[1]}`))));
  t.append(tb);
  $("#products-table").append(el("div", { class: "tablewrap" }, t));
}

function renderBias(bias) {
  const host = $("#bias-table");
  const regions = Object.keys(bias.regions);
  const prods = META.products.map((p) => p.key).filter((k) => k !== bias.ref);
  let variant = "all";
  const draw = () => {
    host.querySelector(".tablewrap")?.remove();
    host.append(heatTable({
      rows: prods, cols: regions, rotate: true,
      value: (r, c) => bias.table[r]?.[c]?.[variant],
      color: ratioColor, format: (v) => fmt(v, 2),
      colLabel: (c) => `${c} (${bias.regions[c].cells})`,
      title: "ratio to GPCC",
    }));
  };
  host.append(el("div", { class: "controls" }, el("label", {}, "Months / footing:"),
    seg([["all", "All months"], ["rain", "Rain months only (T ≥ 2 °C)"], ["lw", "Gauges undercatch-corrected (LW)"]], (v) => { variant = v; draw(); }, "all")),
    scaleLegend(DIV, ["drier than GPCC (ratio < 1)", "wetter"]));
  draw();

  const ph = $("#bias-pairwise");
  const regs = Object.keys(bias.pairwise);
  const keys = META.products.map((p) => p.key);
  let reg = regs[0];
  const drawP = () => {
    ph.querySelector(".tablewrap")?.remove();
    ph.append(heatTable({
      rows: keys, cols: keys, rotate: true,
      value: (r, c) => (r === c ? null : bias.pairwise[reg][r][c]),
      color: ratioColor, format: (v) => fmt(v, 2), colLabel: (c) => SHORT[c],
      title: "row ÷ column",
    }));
  };
  ph.append(el("div", { class: "controls" }, el("label", {}, "Region:"), selectBox(regs.map((r) => [r, r]), (v) => { reg = v; drawP(); }, reg)));
  drawP();

  const uc = bias.undercatch_check;
  if (uc) {
    $("#undercatch-text").innerHTML = `Across the common domain the cell-level CHIRPS v3 / v2 ratio and the annual
      Legates–Willmott factor correlate at <b>r = ${fmt(uc.corr, 2)}</b> (median ratio ${fmt(uc.median_v3_over_v2, 3)}
      vs median factor ${fmt(uc.median_lw, 3)}). Median v3/v2 by factor band: ${uc.binned.map((b) => `${b.lw_lo}–${b.lw_hi}: ${fmt(b.median_v3_over_v2, 3)} (n=${b.n})`).join("; ")}.`;
  }
}

function renderAgreement(ag) {
  const host = $("#agree-table");
  const regions = Object.keys(ag.table[META.products[0].key] || ag.table[Object.keys(ag.table)[0]]);
  const prods = META.products.map((p) => p.key).filter((k) => k !== ag.ref);
  let metric = "spearman_median";
  const draw = () => {
    host.querySelector(".tablewrap")?.remove();
    const isGamma = metric === "variability_ratio_median";
    host.append(heatTable({
      rows: prods, cols: regions, rotate: true,
      value: (r, c) => ag.table[r]?.[c]?.[metric],
      color: (v) => (isGamma ? ratioColor(v) : seqColor(v, 0.2, 1)), format: (v) => fmt(v, 2),
      title: isGamma ? "anomaly s.d. ÷ GPCC's" : "median Spearman ρ",
    }));
  };
  host.append(el("div", { class: "controls" }, el("label", {}, "Metric:"),
    seg([["spearman_median", "Spearman ρ (median cell)"], ["spearman_p25", "Spearman ρ (25th pct cell)"], ["variability_ratio_median", "Variability ratio"]], (v) => { metric = v; draw(); }, metric)));
  draw();
  const ph = $("#agree-pairwise");
  const regs = Object.keys(ag.pairwise_spearman);
  const keys = META.products.map((p) => p.key);
  let reg = regs[0];
  const drawP = () => {
    ph.querySelector(".tablewrap")?.remove();
    ph.append(heatTable({ rows: keys, cols: keys, rotate: true, value: (r, c) => (r === c ? null : ag.pairwise_spearman[reg][r][c]),
      color: (v) => seqColor(v, 0.2, 1), format: (v) => fmt(v, 2), colLabel: (c) => SHORT[c], title: "median ρ" }));
  };
  ph.append(el("div", { class: "controls" }, el("label", {}, "Region:"), selectBox(regs.map((r) => [r, r]), (v) => { reg = v; drawP(); }, reg)));
  drawP();
}

function renderTC(tc) {
  const prods = META.products.map((p) => p.key).filter((k) => tc.table[k]);
  const regions = Object.keys(tc.table[prods[0]].regions);
  $("#tc-table").append(heatTable({
    rows: prods, cols: regions, rotate: true,
    value: (r, c) => tc.table[r].regions[c]?.rho2_median,
    color: (v) => seqColor(v, 0.2, 1), format: (v) => fmt(v, 2),
    rowLabel: (r) => `${LABEL[r]}  [vs ${tc.table[r].partners.map((p) => SHORT[p]).join(", ")}]`,
    title: "median ρ² with truth",
  }));
}

function renderTerciles(tt) {
  const host = $("#terc-table");
  const prods = META.products.map((p) => p.key).filter((k) => k !== tt.ref);
  const regions = Object.keys(tt.table[prods[0]]["3-month"]);
  let scale = "3-month", clim = "own", metric = "pss";
  const draw = () => {
    host.querySelector(".tablewrap")?.remove();
    host.append(heatTable({
      rows: prods, cols: regions, rotate: true,
      value: (r, c) => tt.table[r]?.[scale]?.[c]?.[clim]?.[metric],
      color: (v) => (metric === "pofd" ? seqColor(0.6 - v, 0, 0.6) : seqColor(v, metric === "hit" ? 0.3 : 0, 1)),
      format: (v) => fmt(v, 2), title: { pss: "Peirce skill", kappa: "Cohen's κ (3 classes)", hit: "dry hit rate", pofd: "false-alarm rate" }[metric],
    }));
  };
  host.append(el("div", { class: "controls" },
    el("label", {}, "Accumulation:"), seg([["1-month", "1 month"], ["3-month", "3 months"]], (v) => { scale = v; draw(); }, scale),
    el("label", {}, "Terciles from:"), seg([["own", "each product's own climatology"], ["common", "GPCC's climatology"]], (v) => { clim = v; draw(); }, clim),
    el("label", {}, "Score:"), seg([["pss", "PSS"], ["kappa", "κ"], ["hit", "Hit rate"], ["pofd", "False alarms"]], (v) => { metric = v; draw(); }, metric)));
  draw();
  const ph = $("#terc-pairwise");
  const regs = Object.keys(tt.pairwise_pss_3month);
  const keys = META.products.map((p) => p.key);
  let reg = regs[0];
  const drawP = () => {
    ph.querySelector(".tablewrap")?.remove();
    ph.append(heatTable({ rows: keys, cols: keys, rotate: true, value: (r, c) => (r === c ? null : tt.pairwise_pss_3month[reg][r][c]),
      color: (v) => seqColor(v, 0, 1), format: (v) => fmt(v, 2), colLabel: (c) => SHORT[c], title: "PSS: row detects column's dry terciles" }));
  };
  ph.append(el("div", { class: "controls" }, el("label", {}, "Region:"), selectBox(regs.map((r) => [r, r]), (v) => { reg = v; drawP(); }, reg)));
  drawP();
}

function renderTrends(tr) {
  const host = $("#trend-table");
  let win = "long";
  const draw = () => {
    host.querySelector(".tablewrap")?.remove();
    const W = tr[win];
    const prods = META.products.map((p) => p.key).filter((k) => W.products[k]);
    const regions = [...new Set(prods.flatMap((p) => Object.keys(W.products[p].regions)))];
    host.append(heatTable({
      rows: prods, cols: regions, rotate: true,
      value: (r, c) => W.products[r].regions[c],
      color: (v) => trendColor(v?.pct_dec), format: (v) => (v ? `${pct(v.pct_dec, 1)}${v.p < 0.05 ? " •" : ""}` : "–"),
      title: `% per decade, ${W.years[0]}–${W.years[1]}`,
    }));
    drawDiff();
  };
  host.append(el("div", { class: "controls" }, el("label", {}, "Window:"),
    seg([["long", `${tr.long.years[0]}–${tr.long.years[1]}`], ["recent", `${tr.recent.years[0]}–${tr.recent.years[1]}`]], (v) => { win = v; draw(); }, win)),
    scaleLegend(DIV, ["decreasing", "increasing (• = p < 0.05, modified Mann–Kendall)"]));
  const dh = $("#diff-table");
  function drawDiff() {
    dh.innerHTML = "";
    const D = tr.difference[win];
    const pairs = Object.keys(D);
    const regions = [...new Set(pairs.flatMap((p) => Object.keys(D[p]).filter((k) => !k.startsWith("_"))))];
    dh.append(heatTable({
      rows: pairs, cols: regions, rotate: true,
      value: (r, c) => D[r][c],
      color: (v) => trendColor(v?.pct_dec_of_ref), format: (v) => (v ? `${pct(v.pct_dec_of_ref, 1)}${v.p < 0.05 ? " •" : ""}` : "–"),
      rowLabel: (r) => { const [a, b] = r.split("|"); return `${SHORT[a]} − ${SHORT[b]}`; },
      title: "trend of A − B, % of B per decade",
    }));
  }
  draw();
}

function renderNetworks(nw) {
  const n = nw.networks, yrs = n.years;
  const host = $("#network-chart");
  const S = [];
  const add = (key, label, color, dash, vals, scale = 1) => vals && S.push({ key, label, color, dash, values: vals.map((v) => (v == null ? null : v * scale)) });
  add("gpcc_full_frac_cells_gauged", "GPCC Full: % of cells with ≥1 gauge", "#1f2324", "", n.gpcc_full_frac_cells_gauged, 100);
  add("gpcc_monitoring_frac_cells_gauged", "GPCC Monitoring: % of cells with ≥1 gauge (1° parent)", "#1f2324", "4 3", n.gpcc_monitoring_frac_cells_gauged, 100);
  add("cru_frac_cells_with_station", "CRU TS: % of cells with a station in range", "#9aa3a4", "", n.cru_frac_cells_with_station, 100);
  add("imerg_final_gauge_weight", "IMERG Final: mean gauge weight (%)", "#4a3aa7", "", n.imerg_final_gauge_weight, 1);
  add("chirps_v3_fill_frac", "CHIRPS v3: % of pentads filled with ERA5", "#1baf7a", "5 3", n.chirps_v3_fill_frac, 100);
  lineChart(host, { x: yrs, series: S, yLabel: "%", yFormat: (v) => fmt(v, 0) });
  const h2 = $("#station-chart");
  const S2 = [];
  if (n.gpcc_full_gauges) S2.push({ key: "g", label: "GPCC Full Data: gauges in land cells (monthly mean)", color: "#1f2324", dash: "", values: n.gpcc_full_gauges });
  if (n.chirps_v3_stations) S2.push({ key: "c", label: "CHIRPS v3: stations (monthly mean)", color: "#1baf7a", dash: "", values: n.chirps_v3_stations });
  lineChart(h2, { x: yrs, series: S2, yLabel: "stations", yFormat: (v) => (v >= 1000 ? `${fmt(v / 1000, 0)}k` : fmt(v, 0)) });

  const dh = $("#drift-chart");
  const refs = Object.keys(nw.drift);
  let ref = refs[0], reg = Object.keys(nw.drift[ref])[0];
  const draw = () => {
    const D = nw.drift[ref][reg];
    const keys = META.products.map((p) => p.key).filter((k) => D[k] && k !== ref);
    lineChart($("#drift-chart-plot"), { x: yrs, series: seriesFor(keys, (k) => D[k]), yLabel: `annual total ÷ ${SHORT[ref]}`, yFormat: (v) => fmt(v, 2), refY: 1 });
  };
  dh.append(el("div", { class: "controls" },
    el("label", {}, "Reference:"), selectBox(refs.map((r) => [r, LABEL[r]]), (v) => { ref = v; draw(); }, ref),
    el("label", {}, "Region:"), selectBox(Object.keys(nw.drift[ref]).map((r) => [r, r]), (v) => { reg = v; draw(); }, reg)),
    el("div", { id: "drift-chart-plot" }));
  draw();
  const st = $("#steps-table");
  const t = el("table", { class: "t" });
  t.append(el("thead", {}, el("tr", {}, ...["Test", "Before", "Ratio", "After", "Ratio", "Change"].map((h) => el("th", {}, h)))));
  const tb = el("tbody");
  Object.entries(nw.steps).forEach(([k, v]) => {
    if (!v) return;
    const td = el("td", {}, pct(v.change_pct, 1));
    cellStyle(td, trendColor(v.change_pct));
    tb.append(el("tr", {}, el("td", { class: "rowh" }, k), el("td", {}, v.before.join("–")), el("td", {}, fmt(v.ratio_before, 3)),
      el("td", {}, v.after.join("–")), el("td", {}, fmt(v.ratio_after, 3)), td));
  });
  t.append(tb);
  st.append(el("div", { class: "tablewrap" }, t));
}

function renderCountries(rows) {
  const host = $("#country-table");
  const keys = META.products.map((p) => p.key);
  let metric = "bias";
  const metrics = [["bias", "Bias vs GPCC (ratio)"], ["spearman", "Spearman ρ vs GPCC"], ["pss3", "Dry-tercile PSS (3-month)"],
    ["trend_long", "Trend 1983–2020 (%/dec)"], ["trend_recent", "Trend 2001–2025 (%/dec)"], ["clim", "Mean annual (mm)"]];
  let sortKey = "country", sortDir = 1, filter = "";
  const colorFor = (m, v) => (m === "bias" ? ratioColor(v) : m.startsWith("trend") ? trendColor(v) : m === "clim" ? seqColor(v, 0, 2500) : seqColor(v, m === "pss3" ? 0 : 0.2, 1));
  const fmtFor = (m, v) => (m.startsWith("trend") ? pct(v, 1) : m === "clim" ? fmt(v, 0) : fmt(v, 2));
  const draw = () => {
    host.querySelector(".tablewrap")?.remove();
    const t = el("table", { class: "t sortable" });
    const cols = [["country", "Country"], ["gauges_per_cell", "GPCC gauges / cell"], ["frac_gauge_free", "% cells never gauged"], ...keys.map((k) => [k, SHORT[k]])];
    const hr = el("tr");
    cols.forEach(([k, lab]) => {
      const th = el("th", keys.includes(k) ? { class: "rot" } : {}, keys.includes(k) ? el("div", {}, lab) : lab);
      th.addEventListener("click", () => { sortDir = sortKey === k ? -sortDir : 1; sortKey = k; draw(); });
      hr.append(th);
    });
    t.append(el("thead", {}, hr));
    const val = (r, k) => (k === "country" ? r.country : k === "gauges_per_cell" || k === "frac_gauge_free" ? r[k] : r.products[k]?.[metric]);
    const data = rows.filter((r) => !filter || r.country.toLowerCase().includes(filter) || r.iso3.toLowerCase().includes(filter));
    data.sort((a, b) => {
      const va = val(a, sortKey), vb = val(b, sortKey);
      if (va == null) return 1; if (vb == null) return -1;
      return (typeof va === "string" ? va.localeCompare(vb) : va - vb) * sortDir;
    });
    const tb = el("tbody");
    data.forEach((r) => {
      const tr = el("tr", { style: "cursor:pointer" });
      tr.addEventListener("click", () => explore(r));
      tr.append(el("td", { class: "rowh" }, `${r.country}`), el("td", {}, fmt(r.gauges_per_cell, 1)), el("td", {}, pct(r.frac_gauge_free * 100, 0).replace("+", "")));
      keys.forEach((k) => {
        const v = r.products[k]?.[metric];
        const p = r.products[k]?.[`${metric}_p`];
        const td = el("td", {}, `${fmtFor(metric, v)}${metric.startsWith("trend") && p != null && p < 0.05 ? " •" : ""}`);
        cellStyle(td, colorFor(metric, v));
        tr.append(td);
      });
      tb.append(tr);
    });
    t.append(tb);
    host.append(el("div", { class: "tablewrap", style: "max-height:560px;overflow:auto" }, t));
  };
  const search = el("input", { type: "search", placeholder: "Filter countries…", style: "font:inherit;font-size:13px;padding:4px 8px;border:1px solid #c9cfcf;border-radius:4px" });
  search.addEventListener("input", () => { filter = search.value.toLowerCase(); draw(); });
  host.append(el("div", { class: "controls" }, el("label", {}, "Show:"), selectBox(metrics, (v) => { metric = v; draw(); }, metric), search));
  draw();

  async function explore(r) {
    const s = await getJSON(`country/${r.iso3}`);
    $("#explorer-title").textContent = `${r.country} — ${r.cells} cells, ${fmt(r.gauges_per_cell, 2)} GPCC gauges per cell (mean 2001–2020), ${fmt(r.frac_gauge_free * 100, 0)}% of cells never gauged`;
    const keysA = keys.filter((k) => s.annual[k]);
    lineChart($("#explorer-annual"), { x: s.years, series: seriesFor(keysA, (k) => s.annual[k]), yLabel: "annual total (mm)" });
    lineChart($("#explorer-clim"), { x: ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"], series: seriesFor(keysA, (k) => s.clim[k]), yLabel: "mm/month, 2001–2020", height: 240 });
    $("#explorer").scrollIntoView({ behavior: "smooth", block: "start" });
    history.replaceState(null, "", `#c=${r.iso3}`);
  }
  const m = location.hash.match(/c=([A-Z]{3})/);
  const first = rows.find((r) => r.iso3 === (m ? m[1] : "ETH")) || rows[0];
  if (first) explore(first).then(() => { if (!m) window.scrollTo(0, 0); });
}

// ----------------------------------------------------------------- boot
(async function main() {
  try {
    META = await getJSON("meta");
    META.products.forEach((p) => { LABEL[p.key] = p.label; SHORT[p.key] = p.label.replace(/ \(.*\)$/, "").replace(" Data v2022", "").replace(" v2022", ""); });
    renderProducts(META);
    const [bias, ag, tc, tt, tr, nw, ct] = await Promise.all(["bias", "agreement", "tc", "terciles", "trends", "networks", "countries"].map(getJSON));
    renderBias(bias); renderAgreement(ag); renderTC(tc); renderTerciles(tt); renderTrends(tr); renderNetworks(nw); renderCountries(ct);
  } catch (e) {
    console.error(e);
    const f = $("#fatal"); f.hidden = false; f.textContent = `Could not load report data (${e.message}).`;
  }
})();
