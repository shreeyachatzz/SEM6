highp vec2 func(highp vec2 a, highp vec2 b) {

  // Calculate the length of vectors a and b
  highp float lengthA = length(a);
  highp float lengthB = length(b);
  
  // Calculate the bisector by normalizing the sum of scaled vectors a and b
  return normalize(lengthB * a + lengthA * b);
}

//Do not change this line
#pragma glslify: export(func)