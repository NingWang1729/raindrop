// Raindrop viewer: wires the dataset, the point cloud and the side panel together.
// Open as index.html?data=<dataset folder URL>.
import { Dataset } from './data.js';
import { PointCloudView, SharedCanvas } from './scene.js';
import * as filters from './filters.js';
import { bounds, faceView, toXYZ, validMask } from './positions.js';
import { clampDegrees, faceRotation, openingTurnAxis } from './rotation.js';
import { buildAxisLabels, buildAxisWidget, updateAxisLabels, updateAxisWidget, WIDGET_HALF_SIZE } from './axes.js';
import { Categories } from './categories.js';
import { Coloring } from './coloring.js';
import { Legend } from './legend.js';
import { FilterPanels } from './panels.js';
import { Turning } from './turning.js';
import { el, formatCount, svg } from './ui.js';

// Per-cell states for PointCloudView.setState (see scene.js): not drawn (no position in this embedding), drawn
// faint grey as context (filtered out), drawn in its color (shown)
const NO_POSITION_STATE = 0;
const CONTEXT_STATE = 1;
const SHOWN_STATE = 2;

class App {
  constructor(dataset) {
    this.dataset = dataset;
    this.manifest = dataset.manifest;
    // All views draw into one canvas behind their tiles; see scene.js
    this.canvas = new SharedCanvas(document.getElementById('view'));
    // The views in the grid; the side panel's view and color controls apply to the selected one
    this.views = [];
    this.selected = null;
    // Category colors, the same in every view, and the rows of filter panels and the legend
    this.categories = new Categories(this.manifest.hierarchies);
    // The color menu and gene box, and the colors of each view
    this.coloring = new Coloring({
      manifest: this.manifest,
      dataset,
      categories: this.categories,
      selectedView: () => this.selected,
      setStatus: (message) => this.setStatus(message),
      renderLegend: () => this.legend.renderLegend(),
    });
    // The filter panels and the legend show the selected view's filters; changing one refreshes the views
    this.panels = new FilterPanels({
      manifest: this.manifest,
      dataset,
      categories: this.categories,
      selectedFilters: () => this.settings.filters,
      selectedView: () => this.selected,
      refresh: () => this.refresh(),
    });
    this.legend = new Legend({
      manifest: this.manifest,
      dataset,
      categories: this.categories,
      selectedView: () => this.selected,
      selectedFilters: () => this.settings.filters,
      addFilter: (name) => this.panels.addFilter(name),
      refresh: () => this.refresh(),
    });
    // The side buttons, Turn about menu, rotation table and reset button of the selected view, and the angles shown
    this.turning = new Turning({
      manifest: this.manifest,
      selectedView: () => this.selected,
      selectedSettings: () => this.settings,
      viewsWith: (settings) => this.viewsWith(settings),
      applySettings: (view) => this.applySettings(view),
      showFace: (view, axis, side) => this.showFace(view, axis, side),
      referenceFace: (entry) => this.referenceFace(entry),
    });
  }

  // Build the controls, then open one view on the dataset's default embedding, filters and color
  async init() {
    const { manifest } = this;
    document.title = `${manifest.name} - Raindrop`;
    document.getElementById('title').textContent = manifest.name;
    this.buildViewControls();
    this.coloring.buildColorControls();
    this.legend.buildLegend();
    this.panels.listenForAddFilter();

    const view = this.createView();
    this.select(view);
    await this.setEmbedding(view, manifest.defaults.embedding);
    // Panels for convert.py --filters
    for (const name of manifest.defaults.filters) await this.panels.addFilter(name);
    // Color by convert.py --color, else the first categorical column
    const firstCategory = manifest.fields.find((f) => f.kind === 'category');
    await this.coloring.setColor(view, { kind: 'field', name: manifest.defaults.color ?? firstCategory?.name });
  }

  // The selected view's settings (the shared settings before any view exists)
  get settings() {
    return this.selected?.settings ?? this.group;
  }

  // The views that use a settings object
  viewsWith(settings) {
    return this.views.filter((view) => view.settings === settings);
  }

  // View: embedding, sides, rotation, axes, point size and the opacity of shown and filtered-out cells. Axes, point
  // size, opacities and the turn axis belong to the selected view's settings
  buildViewControls() {
    // Settings shared by the views (filters, turn axis, point size, opacities, axes). Made here, from the
    // controls' starting values
    this.group = {
      filters: [],
      // The axis drags turn about (index into the axes, null: any), from convert.py --locked
      turnAxis: openingTurnAxis(this.manifest.axes.labels, this.manifest.defaults.locked),
      // The point size slider is log2 of the size (0: the default size)
      pointSize: 2 ** Number(document.getElementById('point-size').value),
      shownOpacity: Number(document.getElementById('shown-opacity').value),
      contextOpacity: Number(document.getElementById('context-opacity').value),
      showAxes: document.getElementById('show-axes').checked,
    };

    const select = document.getElementById('embedding');
    for (const e of this.manifest.embeddings) select.append(el('option', { value: e.name, textContent: `${e.name} (${e.dims}D)` }));
    select.addEventListener('change', () => this.setEmbedding(this.selected, select.value));

    this.turning.buildRotationSliders();

    // Each control changes the selected view's settings, and so every view that shares them
    const control = (id, key, read) => {
      const input = document.getElementById(id);
      input.addEventListener(input.type === 'checkbox' ? 'change' : 'input', () => {
        this.settings[key] = read(input);
        for (const view of this.viewsWith(this.settings)) this.applySettings(view);
      });
    };
    control('point-size', 'pointSize', (input) => 2 ** Number(input.value));
    control('shown-opacity', 'shownOpacity', (input) => Number(input.value));
    control('context-opacity', 'contextOpacity', (input) => Number(input.value));
    control('show-axes', 'showAxes', (input) => input.checked);
  }

  // Show a view's display settings and turn axis in it
  applySettings(view) {
    const { settings, scene } = view;
    scene.turnAxis = settings.turnAxis;
    scene.setPointSize(settings.pointSize);
    scene.setShownOpacity(settings.shownOpacity);
    scene.setContextOpacity(settings.contextOpacity);
    scene.setAxesVisible(settings.showAxes);
    view.axisLabelsElement.hidden = !settings.showAxes;
    // With a turn axis, the other axes of a 3D view are drawn as locked (2D views do not turn)
    const turning = view.entry?.dims !== 2 && settings.turnAxis !== null;
    view.widgetAxes.forEach(({ group }, axis) => {
      group.classList.toggle('locked', turning && axis !== settings.turnAxis);
    });
  }

  // A view: its tile (an element over the shared canvas that takes its mouse input, holding the axis labels and
  // the axis widget), its scene, and what it shows (embedding, colors), filled in by setEmbedding and setColor
  createView() {
    const axisLabelsElement = el('div', { className: 'axis-labels' });
    const half = WIDGET_HALF_SIZE;
    const widgetElement = svg('svg', {
      class: 'axis-widget', viewBox: `${-half} ${-half} ${2 * half} ${2 * half}`, 'aria-label': 'Orientation of the data axes',
    });
    const tile = el('div', { className: 'tile' }, [axisLabelsElement, widgetElement]);
    document.getElementById('tiles').append(tile);

    const view = { tile, axisLabelsElement, widgetElement, color: null, settings: this.group };
    view.scene = new PointCloudView(tile, this.canvas);
    view.widgetAxes = buildAxisWidget(widgetElement);
    this.applySettings(view);
    view.scene.onRotate = (degrees) => this.turning.viewTurned(view, degrees);
    // The axis widget and labels follow every drawn frame (turns, zoom, pan, resize)
    view.scene.onRender = () => {
      updateAxisWidget(view);
      updateAxisLabels(view, this.manifest);
    };
    tile.addEventListener('pointerdown', () => this.select(view));
    this.views.push(view);
    // The side panel's Turn about menu shows the selected view's turn axis, which a new view takes from its settings
    this.turning.showTurnAxis();
    return view;
  }

  // Make a view the one the side panel's view and color controls apply to
  select(view) {
    if (this.selected === view) return;
    this.selected = view;
  }

  // Show the side at end `side` (-1 low, +1 high) of `axis`: the camera stays at the reference side (with
  // zoom and pan reset) and the data turns by that side's angles from the reference. faceRotation gives
  // radians; the scene takes degrees, clamped to +-MAX_DEGREES (clampDegrees); the selected view's sliders show the
  // angles (showRotation clamps them too)
  showFace(view, axis, side) {
    const target = faceView(axis, side, view.entry.up, view.entry.roles);
    view.scene.viewFrom(view.reference.direction, view.reference.up);
    const degrees = faceRotation(target, view.reference).map((r) => r * 180 / Math.PI);
    view.scene.setRotation(degrees.map(clampDegrees));
    if (view === this.selected) this.turning.showRotation(degrees);
  }

  // The side the camera looks from, where all angles are 0: the face named by --view if this embedding has
  // it, else the low end of x; 2D data is seen from +z, so x runs right and y up
  referenceFace(entry) {
    if (entry.dims === 2) return { axis: 2, side: 1 };
    const view = this.manifest.defaults.view;
    const axis = entry.ends.findIndex((pair) => pair.includes(view));
    if (axis < 0) return { axis: 0, side: -1 };
    return { axis, side: entry.ends[axis].indexOf(view) === 0 ? -1 : 1 };
  }

  // Show embedding `name` in a view: load its positions (cells without one are left out), frame the camera on
  // the reference side with every angle 0, and recolor and refresh the view
  async setEmbedding(view, name) {
    const entry = this.manifest.embeddings.find((e) => e.name === name);
    view.entry = entry;
    const positions = await this.dataset.embedding(name);
    view.valid = validMask(positions, entry.dims);
    const xyz = toXYZ(positions, entry.dims, view.valid);
    view.box = bounds(xyz, view.valid);
    view.scene.setPositions(xyz, view.box, entry.dims === 2);
    view.scene.setUpAxis(entry.up);
    view.scene.setPointSize(view.settings.pointSize);
    buildAxisLabels(view);
    // A 2D embedding's widget shows its two axes only (SVG elements ignore `hidden`, so a class drops the third),
    // none of them greyed by the turn axis
    view.widgetAxes[2].group.classList.toggle('absent', entry.dims === 2);
    this.applySettings(view);
    // The reference side's view: 2D data keeps +y up whatever its manifest says
    const { axis, side } = this.referenceFace(entry);
    view.reference = faceView(axis, side, entry.dims === 2 ? [0, 1, 0] : entry.up, entry.roles);
    if (view === this.selected) {
      document.getElementById('embedding').value = name;
      this.turning.showViewControls(view);
    }
    this.showFace(view, axis, side);
    if (view.color) await this.coloring.setColor(view, view.color);
    this.refresh();
  }

  // After any filter, embedding or view change: recompute which cells pass in each view (cells without a
  // position differ between embeddings), then update the points, and the count, filter panels, add-filter menu
  // and legend of the selected view
  refresh() {
    const ready = this.views.filter((view) => view.valid);
    if (!this.selected?.valid) return;
    for (const view of ready) {
      view.keep = filters.computeKeep(view.valid, view.settings.filters);
      this.updatePoints(view);
    }
    this.updateCount();
    this.panels.updateFilterPanels();
    this.panels.updateAddFilterOptions();
    this.legend.renderLegend();
  }

  // Shown cells in color, filtered-out cells as context, cells without a position not at all
  updatePoints(view) {
    const state = new Uint8Array(view.keep.length);
    for (let i = 0; i < state.length; i++) state[i] = view.valid[i] ? (view.keep[i] ? SHOWN_STATE : CONTEXT_STATE) : NO_POSITION_STATE;
    view.scene.setState(state);
  }

  // The selected view's count: cells passing its filters, of those with a position in its embedding
  updateCount() {
    const nShown = filters.countKept(this.selected.keep);
    const nValid = filters.countKept(this.selected.valid);
    document.getElementById('count').textContent = `${formatCount(nShown)} of ${formatCount(nValid)} cells shown`;
  }

  setStatus(message) {
    document.getElementById('status').textContent = message;
  }
}

// Open the dataset named by ?data=, or say how to; errors are shown in the status line and rethrown for the
// console
async function main() {
  const params = new URLSearchParams(window.location.search);
  const url = params.get('data');
  if (!url) {
    document.getElementById('status').textContent = 'Add ?data=<dataset folder URL> to the address to open a dataset.';
    return;
  }
  try {
    const app = new App(await Dataset.load(url));
    await app.init();
  } catch (error) {
    document.getElementById('status').textContent = `Error: ${error.message}`;
    throw error;
  }
}

main();
