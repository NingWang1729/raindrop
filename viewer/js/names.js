// Names the viewer shows for what the data holds. Plain functions over manifest entries, so they run in node
// tests too.

// What a view is colored by, as its legends name it: an obs column by its name; a gene with its scale ('log scale'
// when the viewer takes the log, or the source's own, e.g. 'log2'), and for genes of a source other than the first
// (X's measured genes) also the source, e.g. 'Pcp4 (imputed, log scale)'. `sources` is the manifest's gene_sources
export function colorLabel(color, sources) {
  if (color.kind !== 'gene') return color.name;
  const scale = sources.find((s) => s.name === color.source).scale ?? 'log scale';
  return color.source === sources[0].name ? `${color.name} (${scale})` : `${color.name} (${color.source}, ${scale})`;
}
