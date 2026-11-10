// Tests for the viewer's pure-logic modules (data, positions, colors, filters, hierarchy, rotation, names, picking).
// Run with: node tests/test_viewer.mjs
import assert from 'node:assert/strict';
import * as filters from '../viewer/js/filters.js';
import { Dataset, decodeSteps } from '../viewer/js/data.js';
import { categoryPalette, colorCategories, colorValues, colormap, generatedColor, geneColors, parseHex, quantile, sortedFinite, toHex } from '../viewer/js/colors.js';
import { bounds, faceView, framedDiameterPx, niceTicks, readsHighToLow, tickTarget, toXYZ, uprightAngle, validMask } from '../viewer/js/positions.js';
import { groupState, treeRows } from '../viewer/js/hierarchy.js';
import { axisLabels, axisTitle, capitalize, cellValueText, colorLabel, endLabel, widgetEndLabel } from '../viewer/js/names.js';
import { buildScreenGrid, makeProjection, nearestCell, projectShownCells } from '../viewer/js/picking.js';
import { anglesToMatrix, clampDegrees, continuousAngles, dragDirection, faceRotation, knobTurn, openingTurnAxis, recentMovement, wrapTurns } from '../viewer/js/rotation.js';

const tests = {
  'category and range filters combine with AND, within the base mask'() {
    const base = Uint8Array.from([1, 1, 1, 1, 1, 0]);
    const section = { kind: 'category', values: Uint8Array.from([0, 0, 1, 1, 2, 0]), allowed: null };
    const depth = { kind: 'range', values: Float32Array.from([1, 5, 2, NaN, 9, 1]), min: null, max: null };

    // Inactive filters keep every base cell
    assert.deepEqual([...filters.computeKeep(base, [section, depth])], [1, 1, 1, 1, 1, 0]);

    // Sections 0 and 1, depth up to 4: NaN fails an active range
    filters.setCategories(section, 3, [0, 1], true);
    depth.max = 4;
    assert.deepEqual([...filters.computeKeep(base, [section, depth])], [1, 0, 1, 0, 0, 0]);

    // Skipping a filter gives the counts shown in that filter's own panel
    const others = filters.computeKeep(base, [section, depth], section);
    assert.deepEqual([...filters.countCategories(section.values, 3, others)], [1, 1, 0]);
    assert.equal(filters.countKept(others), 2);
  },

  'a reversed category filter keeps the values not ticked, and combines with AND'() {
    // Cells in regions 0 and 1, of subclasses 0, 1 and 2; the last cell has no position
    const base = Uint8Array.from([1, 1, 1, 1, 1, 0]);
    const region = { kind: 'category', values: Uint8Array.from([0, 0, 0, 1, 1, 0]), allowed: null, reversed: false };
    const subclass = { kind: 'category', values: Uint8Array.from([0, 1, 2, 0, 1, 0]), allowed: null, reversed: true };

    // Reversed with nothing ticked restricts nothing
    assert.equal(filters.isActive(subclass), false);
    assert.deepEqual([...filters.computeKeep(base, [subclass])], [1, 1, 1, 1, 1, 0]);

    // Region 0 ticked, subclass 0 ticked and reversed: region 0 cells except subclass 0
    filters.setCategories(region, 2, [0], true);
    filters.setCategories(subclass, 3, [0], true);
    assert.deepEqual([...filters.computeKeep(base, [region, subclass])], [0, 1, 1, 0, 0, 0]);
    // The ticks are not rewritten; unreversed, the same ticks show region 0's subclass 0 cells again
    assert.deepEqual([...subclass.allowed], [1, 0, 0]);
    subclass.reversed = false;
    assert.deepEqual([...filters.computeKeep(base, [region, subclass])], [1, 0, 0, 0, 0, 0]);

    // A panel's counts skip its own filter, reversed or not
    subclass.reversed = true;
    const others = filters.computeKeep(base, [region, subclass], subclass);
    assert.deepEqual([...filters.countCategories(subclass.values, 3, others)], [1, 1, 1]);
  },

  'a reversed range keeps the cells outside it, NaN included'() {
    const base = Uint8Array.from([1, 1, 1, 1, 1]);
    const depth = { kind: 'range', values: Float32Array.from([1, 2, 3, 4, NaN]), min: null, max: null, reversed: true };
    // No bounds: no effect, reversed or not
    assert.deepEqual([...filters.computeKeep(base, [depth])], [1, 1, 1, 1, 1]);
    // The closed range [2, 3] reversed keeps 1, 4 and the NaN cell (which fails the range)
    depth.min = 2;
    depth.max = 3;
    assert.deepEqual([...filters.computeKeep(base, [depth])], [1, 0, 0, 1, 1]);
    depth.reversed = false;
    assert.deepEqual([...filters.computeKeep(base, [depth])], [0, 1, 1, 0, 0]);
  },

  'setting several categories at once; none ticked turns the filter off'() {
    const filter = { kind: 'category', values: Uint8Array.from([0, 1, 2]), allowed: null };
    filters.setCategories(filter, 4, [0, 2, 3], true);
    assert.deepEqual([...filter.allowed], [1, 0, 1, 1]);
    filters.setCategories(filter, 4, [2, 3], false);
    assert.deepEqual([...filter.allowed], [1, 0, 0, 0]);
    filters.setCategories(filter, 4, [0], false);
    assert.equal(filter.allowed, null);
  },

  'tree rows: groups with summed counts, children only when expanded or matching a search'() {
    // Subclasses 1 a, 2 b, 10 c, 3 d in classes X (1 a, 10 c) and Y (2 b, 3 d)
    const names = ['1 a', '2 b', '10 c', '3 d'];
    const counts = [5, 1, 2, 0];
    const parentOf = [0, 1, 0, 1];
    const groups = ['X', 'Y'];
    const brief = (rows) => rows.map((r) => (r.kind === 'group' ? `${groups[r.group]}:${r.count}:${r.open}` : names[r.code]));

    assert.deepEqual(brief(treeRows(names, counts, parentOf, groups, 'name', '', new Set())), ['X:7:false', 'Y:1:false']);
    // Expanded groups list their children in natural order; count order sorts groups and children by count
    assert.deepEqual(brief(treeRows(names, counts, parentOf, groups, 'name', '', new Set([0]))), ['X:7:true', '1 a', '10 c', 'Y:1:false']);
    assert.deepEqual(brief(treeRows(names, counts, parentOf, groups, 'count', '', new Set([1]))), ['X:7:false', 'Y:1:true', '2 b', '3 d']);
    // A search opens the groups with matching children and shows only those; a matching group shows all of its children
    assert.deepEqual(brief(treeRows(names, counts, parentOf, groups, 'name', 'c', new Set())), ['X:7:true', '10 c']);
    assert.deepEqual(brief(treeRows(names, counts, parentOf, groups, 'name', 'y', new Set())), ['Y:1:true', '2 b', '3 d']);
    assert.deepEqual(treeRows(names, counts, parentOf, groups, 'name', '', new Set())[0].codes, [0, 2]);
  },

  'a group is ticked when all its children are, partly when some are'() {
    assert.equal(groupState([0, 2], null), 'none');
    assert.equal(groupState([0, 2], Uint8Array.from([1, 0, 1])), 'all');
    assert.equal(groupState([0, 2], Uint8Array.from([0, 1, 1])), 'some');
    assert.equal(groupState([0, 2], Uint8Array.from([0, 1, 0])), 'none');
  },

  'categories sort in natural name order, or by count with ties by name'() {
    const names = ['10 Oligo', '2 Astro', '1 Endo', 'b', 'A'];
    const counts = [5, 9, 5, 0, 1];
    assert.deepEqual(filters.orderCategories(names, counts, 'name'), [2, 1, 0, 4, 3]);
    assert.deepEqual(filters.orderCategories(names, counts, 'count'), [1, 2, 0, 4, 3]);
  },

  'angles stop at one turn either way, and drags wrap past it'() {
    assert.deepEqual([0, 270, 360, 400, -270, -500].map(clampDegrees), [0, 270, 360, 360, -270, -360]);
    // Dragged angles past one turn (±2π) continue from near 0 (a whole turn less): the same orientation
    const turn = 2 * Math.PI;
    assert.deepEqual(wrapTurns([0, turn, -turn]), [0, turn, -turn]);
    assert.ok(Math.abs(wrapTurns([turn + 0.01, 0, 0])[0] - 0.01) < 1e-12);
    assert.ok(Math.abs(wrapTurns([0, 0, -turn - 0.5])[2] + 0.5) < 1e-12);
  },

  'angles read from an orientation stay continuous with the previous ones'() {
    const near = (actual, expected) => assert.ok(actual.every((v, i) => Math.abs(v - expected[i]) < 1e-9), `${actual} vs ${expected}`);
    // Past 180 degrees: -3.1 continues as 2 pi - 3.1 rather than jumping
    near(continuousAngles([0, 0, -3.1], [0, 0, 3.0]), [0, 0, 2 * Math.PI - 3.1]);
    // Past 90 degrees in y, Euler angles switch to the other solution (a + pi, pi - b, c + pi); it is undone
    near(continuousAngles([Math.PI, Math.PI - 1.6, Math.PI], [0, 1.55, 0]), [0, 1.6, 0]);
  },

  'face views place the camera beyond the chosen end, dorsal up (anterior up when looking down)'() {
    // Mouse CCF: x anterior->posterior, y dorsal->ventral, z left->right, so up (dorsal) is -y
    const up = [0, -1, 0];
    assert.deepEqual(faceView(2, -1, up), { direction: [0, 0, -1], up });
    assert.deepEqual(faceView(0, 1, up), { direction: [1, 0, 0], up });
    assert.deepEqual(faceView(1, -1, up), { direction: [0, -1, 0], up: [-1, 0, 0] });
    assert.deepEqual(faceView(2, 1, [0, 1, 0]), { direction: [0, 0, 1], up: [0, 1, 0] });
  },

  'face views follow the axes\' anatomical roles, whatever their order'() {
    // Screen right of a view is up x back (back = direction); `+ 0` turns -0 into 0 for deepEqual
    const right = ({ direction, up }) => [
      up[1] * direction[2] - up[2] * direction[1], up[2] * direction[0] - up[0] * direction[2], up[0] * direction[1] - up[1] * direction[0],
    ].map((v) => v + 0);
    // Another frame: x left->right, y posterior->anterior, z ventral->dorsal, so dorsal up is +z
    const roles = ['LR', 'PA', 'VD'];
    const up = [0, 0, 1];
    // Dorsal and ventral views: anterior (+y) up; seen from dorsal, the right side (+x) is on the right
    const dorsal = faceView(2, 1, up, roles);
    assert.deepEqual(dorsal, { direction: [0, 0, 1], up: [0, 1, 0] });
    assert.deepEqual(right(dorsal), [1, 0, 0]);
    assert.deepEqual(faceView(2, -1, up, roles).up, [0, 1, 0]);
    // Left view: dorsal up and anterior on the left (posterior, -y, to the right); anterior view: dorsal up
    const left = faceView(0, -1, up, roles);
    assert.deepEqual(left, { direction: [-1, 0, 0], up });
    assert.deepEqual(right(left), [0, -1, 0]);
    assert.deepEqual(faceView(1, 1, up, roles), { direction: [0, 1, 0], up });
    // Mouse CCF (AP DV LR, dorsal -y) and its dorsal-up mirror image (AP VD RL, dorsal +y): roles give the same
    // views as the rule for data without roles
    for (const [frameRoles, frameUp] of [[['AP', 'DV', 'LR'], [0, -1, 0]], [['AP', 'VD', 'RL'], [0, 1, 0]]]) {
      for (let axis = 0; axis < 3; axis++) {
        for (const side of [-1, 1]) assert.deepEqual(faceView(axis, side, frameUp, frameRoles), faceView(axis, side, frameUp));
      }
    }
  },

  'the view opens turning about the axis left free by two locked ones, or about any with none locked'() {
    const labels = ['AP', 'DV', 'LR'];
    assert.equal(openingTurnAxis(labels, []), null);
    assert.equal(openingTurnAxis(labels, ['AP', 'LR']), 1);
    assert.equal(openingTurnAxis(labels, ['LR', 'DV']), 0);
    assert.throws(() => openingTurnAxis(labels, ['AP']), /lock two axes/);
    assert.throws(() => openingTurnAxis(labels, ['AP', 'DV', 'LR']), /lock two axes/);
  },

  'the recent movement of a drag follows its last moves, so turning the drag switches its main direction'() {
    let recent = [0, 0];
    for (let i = 0; i < 10; i++) recent = recentMovement(recent, [3, 0.5]);
    assert.equal(dragDirection([0, 0], recent), 'x');
    // A few pixels across do not switch it; moving that way for a while does
    recent = recentMovement(recent, [0, 3]);
    assert.equal(dragDirection([0, 0], recent), 'x');
    for (let i = 0; i < 8; i++) recent = recentMovement(recent, [0, 3]);
    assert.equal(dragDirection([0, 0], recent), 'y');
  },

  'a Shift-drag keeps the main direction of its first few pixels'() {
    assert.equal(dragDirection([10, 10], [12, 11]), null);
    assert.equal(dragDirection([10, 10], [15, 12]), 'x');
    assert.equal(dragDirection([10, 10], [8, 4]), 'y');
    assert.equal(dragDirection([0, 0], [-4, 4]), 'x');
    assert.equal(dragDirection([0, 0], [1, 2], 2), 'y');
  },

  'circling the center turns a knob: counterclockwise on screen is positive, the short way round'() {
    const near = (actual, expected) => assert.ok(Math.abs(actual - expected) < 1e-9, `${actual} vs ${expected}`);
    const center = [100, 100];
    // From right of the center to above it (y down): a quarter turn counterclockwise
    near(knobTurn([150, 100], [100, 50], center), Math.PI / 2);
    near(knobTurn([100, 50], [150, 100], center), -Math.PI / 2);
    // Crossing the left side (atan2 jumps from pi to -pi) still gives the small step
    near(knobTurn([50, 99], [50, 101], center), 2 * Math.atan2(1, 50));
    near(knobTurn([120, 100], [140, 100], center), 0);
  },

  'an axis gets about one tick per 70 pixels on screen, between 2 and 10'() {
    assert.equal(tickTarget(700), 10);
    assert.equal(tickTarget(2000), 10);
    assert.equal(tickTarget(280), 4);
    assert.equal(tickTarget(40), 2);
    assert.equal(tickTarget(0), 2);
  },

  'ticks are round numbers (1, 2 or 5 times a power of ten) inside the range'() {
    assert.deepEqual(niceTicks(0, 10, 5), [0, 2, 4, 6, 8, 10]);
    assert.deepEqual(niceTicks(-3.7, 4.2, 5), [-2, 0, 2, 4]);
    assert.deepEqual(niceTicks(2.48, 6.04, 5), [3, 4, 5, 6]);
    assert.deepEqual(niceTicks(0.11, 0.19, 5), [0.12, 0.14, 0.16, 0.18]);
    assert.deepEqual(niceTicks(5, 5), []);
    // By default about ten ticks: an AP range of -6.6..6.6 mm gets one every 2 mm, DV (-4..4) every 1 mm
    assert.deepEqual(niceTicks(-6.6, 6.6), [-6, -4, -2, 0, 2, 4, 6]);
    assert.deepEqual(niceTicks(-4, 4), [-4, -3, -2, -1, 0, 1, 2, 3, 4]);
  },

  'side views turn the data so that side faces the reference camera, up kept'() {
    const near = (actual, expected) => assert.ok(actual.every((v, i) => Math.abs(v - expected[i]) < 1e-9), `${actual} vs ${expected}`);
    const apply = (m, v) => m.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
    // Mouse CCF: dorsal up is -y; the reference is the left side (camera beyond low z)
    const up = [0, -1, 0];
    const reference = faceView(2, -1, up);
    near(faceRotation(reference, reference), [0, 0, 0]);
    near(faceRotation(faceView(2, 1, up), reference), [0, Math.PI, 0]);
    for (let axis = 0; axis < 3; axis++) {
      for (const side of [-1, 1]) {
        const target = faceView(axis, side, up);
        const turn = anglesToMatrix(faceRotation(target, reference));
        near(apply(turn, target.direction), reference.direction);
        near(apply(turn, target.up), reference.up);
      }
    }
  },

  'palette uses dataset colors, then generated ones, and grey for NA'() {
    const palette = categoryPalette({ categories: ['a', 'b', 'NA'], colors: ['#ff0000', null, null] });
    assert.deepEqual(palette[0], [255, 0, 0]);
    assert.deepEqual(palette[1], generatedColor(1));
    assert.deepEqual(palette[2], [190, 190, 190]);
    assert.deepEqual(parseHex('#0a0B0c'), [10, 11, 12]);
  },

  'continuous colors run from the low to the high end of the colormap'() {
    const rgb = colorValues(Float32Array.from([0, 10, 20, NaN]), 0, 10);
    assert.deepEqual([...rgb.slice(0, 3)], colormap(0));
    assert.deepEqual([...rgb.slice(3, 6)], colormap(1));
    assert.deepEqual([...rgb.slice(6, 9)], colormap(1));
    assert.deepEqual([...rgb.slice(9, 12)], [190, 190, 190]);
  },

  'gene counts are colored on a log scale, genes already on a log scale linearly'() {
    // Half the top value: past the middle of the colormap for counts (log1p), at the middle for log2 values
    const values = Float32Array.from([0, 5, 10]);
    const middle = (rgb) => [...rgb.slice(3, 6)];
    assert.deepEqual(middle(geneColors(values, 10)), colormap(Math.round(255 * Math.log1p(5) / Math.log1p(10)) / 255));
    assert.deepEqual(middle(geneColors(values, 10, 'log2')), colormap(Math.round(255 * 0.5) / 255));
    assert.deepEqual([...geneColors(values, 10, 'log2').slice(6, 9)], colormap(1));
  },

  'one-byte gene values: codes times the gene\'s step, 255 missing'() {
    const values = decodeSteps(Uint8Array.from([0, 127, 254, 255]), 0.5);
    assert.deepEqual([...values], [0, 63.5, 127, NaN]);
    assert.ok(values instanceof Float32Array);
  },

  async 'a dataset downloads each array once, decoding one-byte genes and refusing unknown genes'() {
    // Stand-in for fetch: serves the bytes of a few files and counts the requests
    const files = {
      'http://x/data/genes/0.bin': Uint8Array.from([0, 2, 255]),
      'http://x/data/obs/cluster.bin': Uint16Array.from([3, 1, 2]),
    };
    const requested = [];
    const realFetch = globalThis.fetch;
    globalThis.fetch = async (url) => {
      requested.push(String(url));
      const bytes = files[String(url)];
      return bytes ? { ok: true, arrayBuffer: async () => bytes.buffer } : { ok: false, status: 404 };
    };
    try {
      const manifest = {
        n_cells: 3,
        embeddings: [{ name: 'umap', file: 'missing.bin' }],
        fields: [{ name: 'cluster', file: 'obs/cluster.bin', dtype: 'uint16' }],
        gene_sources: [{ name: 'counts', genes: ['Gad1'], path: 'genes/{index}.bin', dtype: 'uint8', steps: [0.5] }],
      };
      const dataset = new Dataset(new URL('http://x/data/'), manifest);
      assert.equal(dataset.n, 3);
      // Codes times the gene's step, 255 missing; asking twice (even concurrently) downloads once
      const [first, second] = await Promise.all([dataset.geneValues('counts', 'Gad1'), dataset.geneValues('counts', 'Gad1')]);
      assert.deepEqual([...first], [0, 1, NaN]);
      assert.equal(first, second);
      // Fields come in the dtype the manifest gives
      const cluster = await dataset.fieldValues('cluster');
      assert.ok(cluster instanceof Uint16Array);
      assert.deepEqual([...cluster], [3, 1, 2]);
      assert.deepEqual(requested, ['http://x/data/genes/0.bin', 'http://x/data/obs/cluster.bin']);
      // A gene not in the source throws at once; a failed download rejects with the URL and status
      assert.throws(() => dataset.geneValues('counts', 'Sst'), /Sst is not in gene source counts/);
      await assert.rejects(dataset.embedding('umap'), /Could not load http:\/\/x\/data\/missing.bin \(404\)/);
    } finally {
      globalThis.fetch = realFetch;
    }
  },

  'category codes take their palette color; colors convert to and from hex'() {
    const palette = [[255, 0, 0], [0, 128, 255]];
    assert.deepEqual([...colorCategories(Uint8Array.from([1, 0, 1]), palette)], [0, 128, 255, 255, 0, 0, 0, 128, 255]);
    assert.equal(toHex([0, 128, 255]), '#0080ff');
    assert.deepEqual(parseHex(toHex([7, 200, 99])), [7, 200, 99]);
  },

  'generated colors are distinct, valid RGB'() {
    const colors = Array.from({ length: 50 }, (_, i) => generatedColor(i).join(','));
    assert.equal(new Set(colors).size, 50);
    for (const c of colors) assert.ok(c.split(',').every((v) => v >= 0 && v <= 255));
  },

  'quantiles of the finite (or positive) values'() {
    const values = Float32Array.from([0, 3, 0, 1, 2, 0, NaN]);
    assert.deepEqual(sortedFinite(values), [0, 0, 0, 1, 2, 3]);
    assert.deepEqual(sortedFinite(values, true), [1, 2, 3]);
    assert.equal(quantile(sortedFinite(values), 0.99), 3);
    assert.equal(quantile(sortedFinite(values, true), 0), 1);
    assert.throws(() => quantile([], 0.5), /empty/);
  },

  'axis names stay upright and read in screen order: an edge pointing left or straight down reads high end first'() {
    // Screen directions (x right, y down) from an edge's low end to its high end
    const cases = [
      // [dx, dy, angle of the name, whether it reads from the high end]
      [1, 0, 0, false], // low end left: e.g. Anterior-Posterior with anterior on the left
      [-1, 0, 0, true], // low end right: turned upright, so it reads from the high end
      [0, -1, -90, false], // low end at the bottom: reads bottom to top, low end first
      [0, 1, -90, true], // low end on top (dorsal up for a dorsal-to-ventral axis): reads from the high end
      [1, 1, 45, false],
      [-1, -1, 45, true],
    ];
    for (const [dx, dy, angle, reversed] of cases) {
      assert.ok(Math.abs(uprightAngle(dx, dy) - angle) < 1e-9, `${dx}, ${dy}: ${uprightAngle(dx, dy)}`);
      assert.equal(readsHighToLow(dx, dy), reversed, `${dx}, ${dy}`);
    }
  },

  'point size readout: the diameter in CSS px at the framed distance, doubling with the size'() {
    // diameter = size * extent / 200 (data units) * tileCssHeight / (2 tan(fov / 2)) (CSS px per data unit at
    // distance 1) / fitDistance. With a 90 degree field of view tan(45) = 1: 1 * 200 / 200 * 800 / 2 / 100 = 4 px
    assert.ok(Math.abs(framedDiameterPx(1, 200, 100, 800, 90) - 4) < 1e-9);
    assert.ok(Math.abs(framedDiameterPx(2, 200, 100, 800, 90) - 8) < 1e-9);
    assert.ok(Math.abs(framedDiameterPx(8, 200, 100, 800, 90) - 32) < 1e-9);
    // Another case: 0.5 * 10 / 200 * 600 / (2 tan(30)) / 20 = 0.6495 px; half as tall a tile gives half the size
    const expected = (0.5 * 10 / 200) * (600 / (2 * Math.tan(Math.PI / 6))) / 20;
    assert.ok(Math.abs(framedDiameterPx(0.5, 10, 20, 600, 60) - expected) < 1e-9);
    assert.ok(Math.abs(framedDiameterPx(0.5, 10, 20, 300, 60) - expected / 2) < 1e-9);
  },

  'cells without a position are masked, kept off the GPU and out of the bounds'() {
    const positions = Float32Array.from([0, 0, 4, 2, NaN, 7]);
    const valid = validMask(positions, 2);
    assert.deepEqual([...valid], [1, 1, 0]);
    const xyz = toXYZ(positions, 2, valid);
    assert.deepEqual([...xyz], [0, 0, 0, 4, 2, 0, 0, 0, 0]);
    assert.deepEqual(bounds(xyz, valid), { center: [2, 1, 0], extent: 4, lo: [0, 0, 0], hi: [4, 2, 0] });
  },

  'legends name a column by its name, and a gene with its scale and, outside the first source, its source'() {
    const sources = [{ name: 'measured', genes: ['Pcp4'] }, { name: 'imputed', genes: ['Pcp4'] }, { name: 'smoothed', genes: ['Pcp4'], scale: 'log2' }];
    assert.equal(colorLabel({ kind: 'field', name: 'region' }, sources), 'region');
    // Counts (no scale in the manifest) are colored on the viewer's log scale
    assert.equal(colorLabel({ kind: 'gene', source: 'measured', name: 'Pcp4' }, sources), 'Pcp4 (log scale)');
    assert.equal(colorLabel({ kind: 'gene', source: 'imputed', name: 'Pcp4' }, sources), 'Pcp4 (imputed, log scale)');
    // A --log-scaled source names its own scale
    assert.equal(colorLabel({ kind: 'gene', source: 'smoothed', name: 'Pcp4' }, sources), 'Pcp4 (smoothed, log2)');
    assert.equal(colorLabel({ kind: 'gene', source: 'measured', name: 'Pcp4' }, [{ name: 'measured', scale: 'log2' }]), 'Pcp4 (log2)');
  },

  'axis ends are named by their first letter, generic ends as they are, the widget names generic axes by their last word, and names are capitalized'() {
    assert.equal(endLabel('anterior'), 'A');
    assert.equal(endLabel('-x'), '-x');
    assert.equal(capitalize('ventral'), 'Ventral');
    // The widget labels an anatomical end by its letter, and a generic axis ('-UMAP 1', '+UMAP 1') at its high end only
    assert.equal(widgetEndLabel('dorsal'), 'D');
    assert.equal(widgetEndLabel('+UMAP 1'), '1');
    assert.equal(widgetEndLabel('-UMAP 1'), '');
    assert.equal(widgetEndLabel('+x'), 'x');
  },

  'axis names spell out both ends in the given order, or use the axis label for generic ends'() {
    assert.equal(axisTitle('AP', ['anterior', 'posterior']), 'Anterior\u2013Posterior');
    assert.equal(axisTitle('DV', ['ventral', 'dorsal']), 'Ventral\u2013Dorsal');
    assert.equal(axisTitle('UMAP 1', ['-UMAP 1', '+UMAP 1']), 'UMAP 1');
    // An embedding's own axis labels win over the dataset's (older manifests have none)
    const manifest = { axes: { labels: ['AP', 'DV', 'LR'] } };
    assert.deepEqual(axisLabels(manifest, { labels: ['UMAP 1', 'UMAP 2', 'UMAP 3'] }), ['UMAP 1', 'UMAP 2', 'UMAP 3']);
    assert.deepEqual(axisLabels(manifest, {}), ['AP', 'DV', 'LR']);
  },

  'hover picking projects only shown cells in front of the camera and inside the tile, to tile px with depth w'() {
    // Data -> clip matrix (column-major, as three.js stores it): clip = (x, y, 0, z), so NDC = (x / z, y / z, 0) and
    // the depth w is z. In a 100 x 100 tile NDC x 0.5 is px 75 and NDC y 0.5 is px 25 (screen y points down)
    const matrix = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];
    const xyz = Float32Array.from([
      1, 1, 2, // shown, NDC (0.5, 0.5)
      0, 0, 2, // filtered out (state 1)
      0, 0, 2, // no position (state 0)
      0, 0, -1, // behind the camera (w < 0)
      3, 0, 1, // outside the tile (NDC x 3)
      0, 0, 4, // shown, at the center
    ]);
    const state = Float32Array.from([2, 1, 0, 2, 2, 2]);
    const projected = projectShownCells(makeProjection(6), xyz, state, 2, matrix, 100, 100);
    assert.equal(projected.count, 2);
    assert.deepEqual([...projected.cells.subarray(0, 2)], [0, 5]);
    assert.deepEqual([...projected.x.subarray(0, 2)], [75, 50]);
    assert.deepEqual([...projected.y.subarray(0, 2)], [25, 50]);
    assert.deepEqual([...projected.depth.subarray(0, 2)], [2, 4]);
  },

  'hover picking takes the cell nearest the camera within the radius, ties going to the nearest to the pointer'() {
    // clip = (x, y, 0, z) as above; in a 100 x 100 tile a cell at depth z and px (sx, 50) has x = (sx / 50 - 1) * z
    const matrix = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];
    const at = (sx, z) => [(sx / 50 - 1) * z, 0, z];
    const pick = (cells, states, px, radius = 4) => {
      const xyz = Float32Array.from(cells.flat());
      const projected = projectShownCells(makeProjection(states.length), xyz, Float32Array.from(states), 2, matrix, 100, 100);
      return nearestCell(projected, buildScreenGrid(null, projected, 100, 100, 4), px, 50, radius);
    };
    // Cell 1, 2 px from the pointer but nearer the camera, wins over cell 0 under the pointer
    assert.equal(pick([at(50, 4), at(52, 2)], [2, 2], 50), 1);
    // Farther than the radius from the pointer: none
    assert.equal(pick([at(50, 4), at(52, 2)], [2, 2], 60), -1);
    // Equal depth (as every cell of a 2D embedding): the one nearer the pointer
    assert.equal(pick([at(49, 3), at(53, 3)], [2, 2], 52), 1);
    // A cell in the neighbouring bucket (buckets are 4 px: px 47.9 and 50.5 fall in different ones) is found
    assert.equal(pick([at(47.9, 3)], [2], 50.5), 0);
    // A filtered-out cell nearer the camera is skipped
    assert.equal(pick([at(50, 4), at(50, 2)], [2, 1], 50), 0);
    // No shown cell at all
    assert.equal(pick([at(50, 4)], [1], 50), -1);
    // A cell on the tile's right edge (px 100) is kept in the last column of buckets, and found
    assert.equal(pick([at(100, 3)], [2], 98), 0);
  },

  'hover picking follows a turn of the data, and reuses its buffers for the next picture'() {
    // Two cells on the line of sight through the tile's center, at data z = +1 and -1. The camera looks down -z from
    // z = 5: camera z = data z - 5, and w = -(camera z) = 5 - z, so the cell at z = +1 is nearer (w 4). Turned by
    // 180 degrees about y (x -> -x, z -> -z) the other one is nearer. clip = (x, y, 0, w); column-major entries
    // m[8 + 3] and m[12 + 3] give w = -z + 5, or w = z + 5 after the turn
    const xyz = Float32Array.from([0, 0, 1, 0, 0, -1]);
    const state = Float32Array.from([2, 2]);
    const front = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 5];
    const turned = [-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 5];
    const projection = makeProjection(2);
    projectShownCells(projection, xyz, state, 2, front, 100, 100);
    const grid = buildScreenGrid(null, projection, 100, 100, 4);
    assert.equal(nearestCell(projection, grid, 50, 50, 4), 0);
    projectShownCells(projection, xyz, state, 2, turned, 100, 100);
    const regrid = buildScreenGrid(grid, projection, 100, 100, 4);
    assert.equal(regrid.order, grid.order);
    assert.equal(nearestCell(projection, regrid, 50, 50, 4), 1);
  },

  'hover labels give a category its name, and a gene or number its name and value to 3 significant digits'() {
    const sources = [{ name: 'measured' }, { name: 'smoothed', scale: 'log2' }];
    const region = { kind: 'category', categories: ['r0', 'r1', 'NA'] };
    assert.equal(cellValueText({ kind: 'field', name: 'region' }, region, 1, sources), 'r1');
    assert.equal(cellValueText({ kind: 'field', name: 'region' }, region, 2, sources), 'NA');
    assert.equal(cellValueText({ kind: 'gene', source: 'measured', name: 'Pcp4' }, null, 2.3456, sources), 'Pcp4: 2.35');
    assert.equal(cellValueText({ kind: 'gene', source: 'measured', name: 'Pcp4' }, null, 7, sources), 'Pcp4: 7');
    assert.equal(cellValueText({ kind: 'field', name: 'depth' }, { kind: 'numeric' }, 1234.5, sources), 'depth: 1230');
    assert.equal(cellValueText({ kind: 'gene', source: 'measured', name: 'Pcp4' }, null, NaN, sources), 'Pcp4: no value');
    // A source on its own scale (--log-scaled) names it, as the value is not a count
    assert.equal(cellValueText({ kind: 'gene', source: 'smoothed', name: 'Pcp4' }, null, 3.14159, sources), 'Pcp4 (log2): 3.14');
  },
};

let failed = 0;
for (const [name, test] of Object.entries(tests)) {
  try {
    await test();
    console.log(`ok   ${name}`);
  } catch (error) {
    failed++;
    console.log(`FAIL ${name}\n${error.stack}`);
  }
}
process.exit(failed ? 1 : 0);
