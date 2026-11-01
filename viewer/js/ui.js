// Small DOM helpers shared by the viewer's modules. Elements are built with properties (textContent, not
// innerHTML), so names that come from the data are never parsed as HTML.

// Create an element with properties and children
export function el(tag, props = {}, children = []) {
  const element = document.createElement(tag);
  Object.assign(element, props);
  element.append(...children);
  return element;
}

// A count with thousands separators, the same in every locale (e.g. 12,345)
export function formatCount(n) {
  return n.toLocaleString('en-US');
}
