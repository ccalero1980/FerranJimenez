
  
float dotProduct(vec_s a, vec_s b)
{
	return a.x*b.x+a.y*b.y+a.z*b.z;
}
 
vec_s crossProduct(vec_s a, vec_s b)
{
	vec_s c = {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
 
	return c;
}

vec_s scalarProduct(vec_s a, float b)
{
	vec_s c = {a.x*b, a.y*b, a.z*b};
 
	return c;
}
 
vec_s vecAddition(vec_s a, vec_s b)
{
	vec_s c = {a.x + b.x, a.y + b.y, a.z + b.z};
 
	return c;
}
 
 
