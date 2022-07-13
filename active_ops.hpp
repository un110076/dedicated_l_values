#pragma once

active operator+(const active& arg1, const active& arg2) {
  active res;
  if (arg1.node_ref>=0&&arg2.node_ref>=0) {
    if (arg1.node_ref==arg2.node_ref) {
      arg1.record_arg(2);
      res.record_res(1);
    } else {
      arg1.record_arg(1.);
      arg2.record_arg(1.);
      res.record_res(2);
    }
  } else if (arg1.node_ref>=0) {
    arg1.record_arg(1.);
    res.record_res(1);
  } else if (arg2.node_ref>=0) {
    arg2.record_arg(1.);
    res.record_res(1);
  } 
  res.value=arg1.value+arg2.value;
  return res;
}

active operator-(const active& arg1, const active& arg2) {
  active res;
  if (arg1.node_ref>=0&&arg2.node_ref>=0) {
    if (arg1.node_ref!=arg2.node_ref) {
      arg1.record_arg(1.);
      arg2.record_arg(-1.);
      res.record_res(2);
    }
  } else if (arg1.node_ref>=0) {
    arg1.record_arg(1.);
    res.record_res(1);
  } else if (arg2.node_ref>=0) {
    arg2.record_arg(-1.);
    res.record_res(1);
  } 
  res.value=arg1.value-arg2.value;
  return res;
}

active operator-(const active& arg) {
  active res;
  if (arg.node_ref>=0) {
    arg.record_arg(-1.);
    res.record_res(1);
  } 
  res.value=-arg.value;
  return res;
}

active operator*(const active& arg1, const active& arg2) {
  active res;
  if (arg1.node_ref>=0&&arg2.node_ref>=0) {
    if (arg1.node_ref==arg2.node_ref) {
      arg1.record_arg(2*arg1.value);
      res.record_res(1);
    } else {
      arg1.record_arg(arg2.value);
      arg2.record_arg(arg1.value);
      res.record_res(2);
    }
  } else if (arg1.node_ref>=0) {
    arg1.record_arg(arg2.value);
    res.record_res(1);
  } else if (arg2.node_ref>=0) {
    arg2.record_arg(arg1.value);
    res.record_res(1);
  } 
  res.value=arg1.value*arg2.value;
  return res;
}

active operator/(const active& arg1, const active& arg2) {
  active res;
  if (arg1.node_ref>=0&&arg2.node_ref>=0) {
    if (arg1.node_ref!=arg2.node_ref) {
      arg1.record_arg(1./arg2.value);
      arg2.record_arg(-arg1.value/pow(arg2.value,2));
      res.record_res(2);
    }
  } else if (arg1.node_ref>=0) {
    arg1.record_arg(1./arg2.value);
    res.record_res(1);
  } else if (arg2.node_ref>=0) {
    arg2.record_arg(-arg1.value/pow(arg2.value,2));
    res.record_res(1);
  } 
  res.value=arg1.value/arg2.value;
  return res;
}

bool operator>(const active& arg1, const active& arg2) {
  return arg1.value>arg2.value;
}

bool operator<(const active& arg1, const active& arg2) {
  return arg1.value<arg2.value;
}

active max(const active& arg1, const active& arg2) {
  assert(std::fabs(arg1.value-arg2.value)>sqrt(std::numeric_limits<double>::epsilon()));
  active res;
  if (arg1.node_ref>=0&&arg2.node_ref>=0) {
    if (arg1.value>arg2.value) 
      arg1.record_arg(1.); 
     else 
      arg2.record_arg(1.);
    res.record_res(1);
  } else if (arg1.node_ref>=0&&arg1.value>arg2.value) {
    arg1.record_arg(1.); 
    res.record_res(1);
  } else if (arg2.node_ref>=0&&arg1.value<arg2.value) {
    arg2.record_arg(1.);
    res.record_res(1);
  } 
  res.value=std::max(arg1.value,arg2.value);
  return res;
}

active fabs(const active& arg) {
  active res;
  res.value=fabs(arg.value);
  if (arg.node_ref>=0) {
    if (arg.value<0) 
      arg.record_arg(-1);
    else
      arg.record_arg(1);
    res.record_res(1);
  }
  return res;
}

active exp(const active& arg) {
  active res;
  res.value=exp(arg.value);
  if (arg.node_ref>=0) {
    arg.record_arg(res.value);
    res.record_res(1);
  }
  return res;
}

active sin(const active& arg) {
  active res;
  if (arg.node_ref>=0) {
    arg.record_arg(cos(arg.value));
    res.record_res(1);
  }
  res.value=sin(arg.value);
  return res;
}

active cos(const active& arg) {
  active res;
  if (arg.node_ref>=0) {
    arg.record_arg(-sin(arg.value));
    res.record_res(1);
  }
  res.value=cos(arg.value);
  return res;
}

active sqrt(const active& arg) {
  active res;
  res.value=sqrt(arg.value);
  if (arg.node_ref>=0) {
    arg.record_arg(1/(2*res.value));
    res.record_res(1);
  }
  return res;
}

active erf(const active& arg) {
  active res;
  if (arg.node_ref>=0) {
    arg.record_arg(2/sqrt(4*atan(1))*exp(-pow(arg.value,2)));
    res.record_res(1);
  }
  res.value=erf(arg.value);
  return res;
}

active log(const active& arg) {
  active res;
  if (arg.node_ref>=0) {
    arg.record_arg(1/arg.value);
    res.record_res(1);
  }
  res.value=log(arg.value);
  return res;
}

active pow(const active& arg1, double arg2) {
  active res;
  if (arg1.node_ref>=0) {
    arg1.record_arg(arg2*pow(arg1.value,arg2-1));
    res.record_res(1);
  }
  res.value=pow(arg1.value,arg2);
  return res;
}
