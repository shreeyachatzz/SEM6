bool inBox(highp vec2 lo, highp vec2 hi, highp vec2 p) {

  // Test if the point p is inside the box bounded by [lo, hi]
  bool withinXBounds = p.x >= lo.x && p.x <= hi.x;
  bool withinYBounds = p.y >= lo.y && p.y <= hi.y;

  return withinXBounds && withinYBounds;
}


//Do not change this line or the name of the above function
#pragma glslify: export(inBox)
