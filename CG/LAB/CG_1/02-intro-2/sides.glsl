void sideLengths(
  highp float hypotenuse, 
  highp float angleInDegrees, 
  out highp float opposite, 
  out highp float adjacent) {


  // Convert angle from degrees to radians
  highp float angleInRadians = radians(angleInDegrees);
  
  // Calculate the side lengths using trigonometric functions
  opposite = hypotenuse * sin(angleInRadians);
  adjacent = hypotenuse * cos(angleInRadians);

}

//Do not change this line
#pragma glslify: export(sideLengths)