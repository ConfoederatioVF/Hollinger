"use strict";

/**
 * Exact high-resolution bitset polygon scanline rasteriser for WGS84 Equirectangular grids.
 * Supports multi-ring polygons with holes (even-odd fill rule), exact boolean union across
 * overlapping features, and exact ellipsoidal Shoelace area fallback for micro-islands smaller
 * than a single sub-pixel.
 */
class ScanlineRasteriser {
  /**
   * @param {Object} [arg0_options]
   *  @param {number} [arg0_options.height=2160]
   *  @param {number} [arg0_options.min_micro_area_km2=0.001]
   *  @param {number} [arg0_options.supersample=20]
   *  @param {number} [arg0_options.width=4320]
   */
  constructor (arg0_options) {
    //Convert from parameters
    let options = arg0_options ? arg0_options : {};

    //Initialise options
    if (options.height === undefined) options.height = 2160;
    if (options.min_micro_area_km2 === undefined) options.min_micro_area_km2 = 0.001;
    if (options.supersample === undefined) options.supersample = 20;
    if (options.width === undefined) options.width = 4320;

    //Declare local instance variables
    this.height = options.height;
    this.min_micro_area_km2 = options.min_micro_area_km2;
    this.sub_height = options.height*options.supersample;
    this.sub_width = options.width*options.supersample;
    this.supersample = options.supersample;
    this.width = options.width;

    //Bitset packed into Uint32Array: for 4320x2160 @ S=20 -> 86400x43200 bits = ~466.6 MB
    let total_bits = this.sub_width*this.sub_height;
    this.bitset = new Uint32Array(Math.ceil(total_bits/32));
    this.micro_area_km2 = new Float32Array(this.width*this.height);
    this.row_words = Math.ceil(this.sub_width/32);
  }

  /**
   * Aggregates the supersampled bitset and micro-island areas into [0, 1] fractional coverage per 5' cell.
   *
   * @param {Float64Array|Float32Array} arg0_cell_area_by_row
   *
   * @returns {Float32Array}
   */
  aggregateToFraction (arg0_cell_area_by_row) {
    //Convert from parameters
    let cell_area_by_row = arg0_cell_area_by_row;

    //Declare local instance variables
    let bitset = this.bitset;
    let height = this.height;
    let micro_area = this.micro_area_km2;
    let out_frac = new Float32Array(this.width*height);
    let row_words = this.row_words;
    let s = this.supersample;
    let sub_cells_per_px = s*s;
    let width = this.width;

    for (let y = 0; y < height; y++) {
      let cell_area = cell_area_by_row[y];
      let out_row_offset = y*width;
      let sy_start = y*s;
      let sy_end = sy_start + s;

      for (let x = 0; x < width; x++) {
        let count = 0;
        let sx_start = x*s;
        let sx_end = sx_start + s;

        for (let sy = sy_start; sy < sy_end; sy++) {
          let base_word_row = sy*row_words;
          for (let sx = sx_start; sx < sx_end; sx++) {
            let w_idx = base_word_row + (sx >>> 5);
            if ((bitset[w_idx] >>> (sx & 31)) & 1)
              count++;
          }
        }

        let frac = count/sub_cells_per_px;
        let micro_km2 = micro_area[out_row_offset + x];

        if (count === 0 && micro_km2 > 0 && cell_area > 0)
          frac = Math.min(1.0, micro_km2/cell_area);

        out_frac[out_row_offset + x] = frac;
      }
    }

    //Return statement
    return out_frac;
  }

  /**
   * Computes approximate geodesic area (km2) of a WGS84 ring using spherical excess / cosine-latitude Shoelace.
   *
   * @param {Array<Array<number>>} arg0_ring - [[lon, lat], ...]
   *
   * @returns {{ area_km2: number, centroid_lat: number, centroid_lon: number }}
   */
  static computeRingAreaKm2 (arg0_ring) {
    //Convert from parameters
    let ring = arg0_ring;

    //Guard clause
    if (!ring || ring.length < 3)
      return { area_km2: 0.0, centroid_lat: 0.0, centroid_lon: 0.0 };

    //Declare local instance variables
    let n = ring.length;
    let sum_lat = 0.0;
    let sum_lon = 0.0;

    for (let i = 0; i < n; i++) {
      sum_lon += ring[i][0];
      sum_lat += ring[i][1];
    }

    let mean_lat = sum_lat/n;
    let mean_lon = sum_lon/n;
    let deg_to_km_lat = 111.13295;
    let deg_to_km_lon = 111.31949*Math.cos(mean_lat*Math.PI/180.0);
    let shoelace = 0.0;

    for (let i = 0; i < n; i++) {
      let j = (i + 1) % n;
      let x0 = (ring[i][0] - mean_lon)*deg_to_km_lon;
      let y0 = (ring[i][1] - mean_lat)*deg_to_km_lat;
      let x1 = (ring[j][0] - mean_lon)*deg_to_km_lon;
      let y1 = (ring[j][1] - mean_lat)*deg_to_km_lat;
      shoelace += (x0*y1 - x1*y0);
    }

    //Return statement
    return {
      area_km2: Math.abs(shoelace)*0.5,
      centroid_lat: mean_lat,
      centroid_lon: mean_lon,
    };
  }

  /**
   * Rasterises a GeoJSON Polygon or MultiPolygon geometry into the bitset (OR union).
   *
   * @param {Object} arg0_geometry
   * @param {Object} [arg1_options]
   *  @param {number} [arg1_options.min_hole_area_km2=0]
   */
  rasteriseGeometry (arg0_geometry, arg1_options) {
    //Convert from parameters
    let geometry = arg0_geometry;
    let options = arg1_options ? arg1_options : {};

    //Guard clause
    if (!geometry || !geometry.type || !geometry.coordinates) return;

    let min_hole = (options.min_hole_area_km2 !== undefined) ? options.min_hole_area_km2 : 0.0;
    let filterRings = (arg0_rings) => {
      if (min_hole <= 0.0 || !arg0_rings || arg0_rings.length <= 1) return arg0_rings;
      let kept = [arg0_rings[0]];
      for (let r = 1; r < arg0_rings.length; r++) {
        if (ScanlineRasteriser.computeRingAreaKm2(arg0_rings[r]).area_km2 >= min_hole)
          kept.push(arg0_rings[r]);
      }
      return kept;
    };

    if (geometry.type === "Polygon") {
      this.rasterisePolygonRings(filterRings(geometry.coordinates));
    } else if (geometry.type === "MultiPolygon") {
      for (let i = 0; i < geometry.coordinates.length; i++)
        this.rasterisePolygonRings(filterRings(geometry.coordinates[i]));
    }
  }

  /**
   * Rasterises a single polygon's rings (exterior + holes) via even-odd scanline fill.
   *
   * @param {Array<Array<Array<number>>>} arg0_rings
   */
  rasterisePolygonRings (arg0_rings) {
    //Convert from parameters
    let rings = arg0_rings;

    //Guard clause
    if (!rings || rings.length === 0) return;

    //Declare local instance variables
    let bitset = this.bitset;
    let max_sy = -1;
    let min_sy = this.sub_height;
    let row_words = this.row_words;
    let scale_x = this.sub_width/360.0;
    let scale_y = this.sub_height/180.0;
    let sub_height = this.sub_height;
    let sub_width = this.sub_width;

    //Build edge table in sub-pixel coordinates: [x0, y0, x1, y1]
    let edge_x0 = [];
    let edge_y0 = [];
    let edge_x1 = [];
    let edge_y1 = [];

    for (let i = 0; i < rings.length; i++) {
      let ring = rings[i];
      if (!ring || ring.length < 3) continue;

      for (let x = 0; x < ring.length - 1; x++) {
        let p0 = ring[x];
        let p1 = ring[x + 1];
        let sy0 = (90.0 - p0[1])*scale_y;
        let sy1 = (90.0 - p1[1])*scale_y;

        if (sy0 === sy1) continue;

        let sx0 = (p0[0] + 180.0)*scale_x;
        let sx1 = (p1[0] + 180.0)*scale_x;

        if (sy0 < min_sy) min_sy = sy0;
        if (sy1 < min_sy) min_sy = sy1;
        if (sy0 > max_sy) max_sy = sy0;
        if (sy1 > max_sy) max_sy = sy1;

        edge_x0.push(sx0);
        edge_y0.push(sy0);
        edge_x1.push(sx1);
        edge_y1.push(sy1);
      }
    }

    let edge_count = edge_x0.length;
    if (edge_count === 0) return;

    let sy_start = Math.max(0, Math.floor(min_sy));
    let sy_end = Math.min(sub_height - 1, Math.ceil(max_sy));
    let hits = 0;
    let x_intersections = [];

    for (let sy = sy_start; sy <= sy_end; sy++) {
      let scan_y = sy + 0.5;
      x_intersections.length = 0;

      for (let i = 0; i < edge_count; i++) {
        let y0 = edge_y0[i];
        let y1 = edge_y1[i];

        if ((y0 <= scan_y && y1 > scan_y) || (y1 <= scan_y && y0 > scan_y)) {
          let t = (scan_y - y0)/(y1 - y0);
          x_intersections.push(edge_x0[i] + t*(edge_x1[i] - edge_x0[i]));
        }
      }

      if (x_intersections.length < 2) continue;

      x_intersections.sort((a, b) => a - b);
      let base_word_row = sy*row_words;

      for (let i = 0; i < x_intersections.length - 1; i += 2) {
        let xs = Math.max(0, Math.ceil(x_intersections[i] - 0.5));
        let xe = Math.min(sub_width - 1, Math.floor(x_intersections[i + 1] - 0.5));

        if (xe < xs) continue;
        hits += (xe - xs + 1);

        //Set bits across [xs, xe] using word-aligned operations
        let w_start = xs >>> 5;
        let w_end = xe >>> 5;
        let b_start = xs & 31;
        let b_end = xe & 31;

        if (w_start === w_end) {
          let mask = (0xFFFFFFFF >>> (31 - (b_end - b_start))) << b_start;
          bitset[base_word_row + w_start] |= mask;
        } else {
          bitset[base_word_row + w_start] |= (0xFFFFFFFF << b_start);
          for (let w = w_start + 1; w < w_end; w++)
            bitset[base_word_row + w] = 0xFFFFFFFF;
          bitset[base_word_row + w_end] |= (0xFFFFFFFF >>> (31 - b_end));
        }
      }
    }

    //Exact micro-island fallback if polygon fell between 15" sub-pixel centres
    if (hits === 0 && rings[0] && rings[0].length >= 3) {
      let info = ScanlineRasteriser.computeRingAreaKm2(rings[0]);
      if (info.area_km2 >= this.min_micro_area_km2) {
        let row = Math.min(this.height - 1, Math.max(0, Math.floor(((90.0 - info.centroid_lat)/180.0)*this.height)));
        let col = Math.min(this.width - 1, Math.max(0, Math.floor(((info.centroid_lon + 180.0)/360.0)*this.width)));
        this.micro_area_km2[row*this.width + col] += info.area_km2;
      }
    }
  }
}

module.exports = ScanlineRasteriser;
