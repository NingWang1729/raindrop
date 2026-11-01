// Colors of categorical obs columns, made once per column and shared by every view, so a category has the same
// color in all of them.
import { categoryPalette } from './colors.js';

export class Categories {
  constructor() {
    // RGB per category code, by column name
    this.palettes = new Map();
  }

  // RGB of each category of a manifest field (colors.categoryPalette: the dataset's own colors where given,
  // generated ones elsewhere), made on first use
  palette(field) {
    if (!this.palettes.has(field.name)) this.palettes.set(field.name, categoryPalette(field));
    return this.palettes.get(field.name);
  }
}
