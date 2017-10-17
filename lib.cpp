#include <map>
#include <string>
#include <iostream>
using namespace std;
#include <cmath>

typedef int VARID;
struct VAR{
	string vname;
	VARID vid;
	static int vid_acc;
	VAR(){
		vid=-1;
		vname="?";
	};
	VAR(string vname):vname(vname){
		vid=vid_acc++;
	};
	void clone(VAR *src) const{
		src->vname=vname;
		src->vid=vid;
	};
	bool eq(VAR &rhs) const{
		return (vid==rhs.vid);
	};
	bool operator<(VAR rhs) const{
		return (vid<rhs.vid);
	}
};
int VAR::vid_acc=0;

#define CONSTANT(x) (new _CONSTANT((x)))
#define IDENTITY(x) (new _IDENTITY((x)))

#define ADD(lhs,rhs) (new _ADD((lhs),(rhs)))
#define SUB(lhs,rhs) (new _SUB((lhs),(rhs)))
#define MUL(lhs,rhs) (new _MUL((lhs),(rhs)))
#define DIV(lhs,rhs) (new _DIV((lhs),(rhs)))
#define POW(lhs,rhs) (new _POW((lhs),(rhs)))
#define LOG(lhs,rhs) (new _LOG((lhs),(rhs)))
#define INV(lhs,rhs) (new _INV((lhs),(rhs)))
#define EXP(lhs,rhs) (new _EXP((lhs),(rhs)))

#define LN(x) (new _LN((x)))

#define COS(x) (new _COS((x)))
#define SIN(x) (new _SIN((x)))
#define TAN(x) (new _TAN((x)))
#define COSH(x) (new _COSH((x)))
#define SINH(x) (new _SINH((x)))
#define TANH(x) (new _TANH((x)))

class _CONSTANT;
class der_f{ public:
	virtual operator string() const =0;
	virtual double value_at(map<VARID, double>) const =0;
	double value_at_origin() const{
		map<VARID, double> empty_vmap;
		return value_at(empty_vmap);
	}
	virtual der_f* _der(VAR) const =0;
	der_f* der(VAR id);
	virtual der_f* clone() const =0;
	virtual ~der_f(){};
	virtual bool has_var(VAR) const =0;
	virtual bool is_constant() const {return false;}//Not so precise
};

class _CONSTANT : public der_f { public:
	double v;
	_CONSTANT(double v):v(v){};
	operator string() const {
		return std::to_string(v);
	}
	der_f* _der(VAR id) const{
		return new _CONSTANT(0);
	};
	der_f* clone() const{
		return new _CONSTANT(v);
	};
	bool has_var(VAR id) const{return false;};
	double value_at(map<VARID, double> vmap) const{
		return v;
	}
	bool is_constant() const{
		return true;
	}
};
der_f* der_f::der(VAR id){
	if(is_constant()||!has_var(id)){
		return new _CONSTANT(0);
	}
	return _der(id);
}

class _IDENTITY : public der_f { public:
	VAR id;
	_IDENTITY(VAR id){
		id.clone(&this->id);
	};
	operator string() const {
		return id.vname;
	}
	der_f* _der(VAR _id) const{
		if(id.eq(_id)){
			return new _CONSTANT(1);
		} else {
			return new _CONSTANT(0);
		}
	};
	der_f* clone() const{
		return new _IDENTITY(id);
	};
	
	bool has_var(VAR _id) const{return id.eq(_id);};
	double value_at(map<VARID, double> vmap) const{
		if(vmap.empty()){
			return 0;
		}
		return vmap[id.vid];
	};
};
class unary_opr : virtual public der_f {public:
	virtual string op_sym() const =0;
	operator string() const{
		return op_sym()+"("+(string)(*ihs)+")";
	}
	der_f *ihs;
	unary_opr(der_f *ihs){
		if(ihs->is_constant()){
			this->ihs=new _CONSTANT(ihs->value_at_origin());
			delete ihs;
		}else{
			this->ihs=ihs;
		}
	}
	bool has_var(VAR id) const{return ihs->has_var(id);};
	virtual ~unary_opr(){
		delete ihs;
	}
	bool is_constant() const{return ihs->is_constant();}
};

class binary_opr : virtual public der_f { public:
	virtual string op_sym() const =0;
	operator string() const {
		return "("+(string)(*lhs)+")"+op_sym()+"("+(string)(*rhs)+")";
	}
	der_f *lhs, *rhs;
	binary_opr(der_f* lhs, der_f* rhs){
		if(lhs->is_constant()){
			this->lhs=new _CONSTANT(lhs->value_at_origin());
			delete lhs;
		}else{
			this->lhs=lhs;
		}
		if(rhs->is_constant()){
			this->rhs=new _CONSTANT(rhs->value_at_origin());
			delete rhs;
		}else{
			this->rhs=rhs;
		}
	};
	virtual ~binary_opr(){
		delete lhs;
		delete rhs;
	};
	bool has_var(VAR id) const{
		return (lhs->has_var(id))||(rhs->has_var(id));
	};
	bool is_constant() const{return (lhs->is_constant())&&(rhs->is_constant());}
};
class _ADD : public binary_opr{ public:
	string op_sym() const{return "+";}
	_ADD(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* _der(VAR id) const{
		return new _ADD(lhs->der(id),rhs->der(id));
	};
	der_f* clone() const{
		return new _ADD(lhs->clone(),rhs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return (lhs->value_at(vmap))+(rhs->value_at(vmap));
	}
};
class _SUB : public binary_opr{ public:
	string op_sym() const{return "-";}
	_SUB(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* _der(VAR id) const{
		return new _SUB(lhs->der(id),rhs->der(id));
	};
	der_f* clone() const{
		return new _SUB(lhs->clone(),rhs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return (lhs->value_at(vmap))-(rhs->value_at(vmap));
	}
};
class _MUL : public binary_opr{ public:
	string op_sym() const{return "*";}
	_MUL(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* _der(VAR id) const{
		return new _ADD(
			new _MUL(lhs->der(id),rhs->clone()),
			new _MUL(lhs->clone(),rhs->der(id))
		);
	};
	der_f* clone() const{
		return new _MUL(lhs->clone(),rhs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return (lhs->value_at(vmap))*(rhs->value_at(vmap));
	}
};
class _DIV : public binary_opr{ public:
	string op_sym() const{return "/";}
	_DIV(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* _der(VAR id) const{
		return new _DIV(
			new _SUB(
				new _MUL(lhs->der(id),rhs->clone()),
				new _MUL(lhs->clone(),rhs->der(id))
			),
			new _MUL(
				rhs->clone(),
				rhs->clone()
			)
		);
	};
	der_f* clone() const{
		return new _DIV(lhs->clone(),rhs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return (lhs->value_at(vmap))/(rhs->value_at(vmap));
	}
};

class _NEG: public unary_opr {public:
	string op_sym() const {return "-";}
	_NEG(der_f *ihs):unary_opr(ihs){}
	der_f* _der(VAR id) const{
		return new _NEG(ihs->der(id));
	}
	der_f* clone() const{
		return new _NEG(ihs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return -(ihs->value_at(vmap));
	}
};
class _INV: public unary_opr {public:
	string op_sym() const {return "inv";}
	_INV(der_f *ihs):unary_opr(ihs){}
	der_f* _der(VAR id) const{
		return new _NEG(
			new _DIV(
				ihs->der(id),
				new _MUL(
					ihs->clone(),
					ihs->clone()
				)
			)
		);
	}
	der_f* clone() const{
		return new _INV(ihs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return 1.0/(ihs->value_at(vmap));
	}
};
class _EXP: public unary_opr {public:
	string op_sym() const {return "exp";}
	_EXP(der_f *ihs):unary_opr(ihs){}
	der_f* _der(VAR id) const{
		return new _MUL(
			this->clone(),
			ihs->der(id)
		);
	}
	der_f* clone() const{
		return new _EXP(ihs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return exp(ihs->value_at(vmap));
	};
	~_EXP(){};
};
class _COS;
class _SIN;
class _COS: public unary_opr {public:
	string op_sym() const {return "cos";}
	_COS(der_f *ihs):unary_opr(ihs){}
	der_f* _der(VAR id) const;
	der_f* clone() const{
		return new _COS(ihs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return cos(ihs->value_at(vmap));
	};
	~_COS(){};
};
class _SIN: public unary_opr {public:
	string op_sym() const {return "sin";}
	_SIN(der_f *ihs):unary_opr(ihs){}
	der_f* _der(VAR id) const;
	der_f* clone() const{
		return new _SIN(ihs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return sin(ihs->value_at(vmap));
	};
	~_SIN(){};
};
der_f* _COS::_der(VAR id) const{
	return new _MUL(
		new _NEG(
			new _SIN(
				ihs->clone()
			)
		),
		ihs->der(id)
	);
};
der_f* _SIN::_der(VAR id) const{
	return new _MUL(
		new _COS(ihs->clone()),
		ihs->der(id)
	);
};
class _TAN: public unary_opr {public:
	string op_sym() const {return "tan";}
	_TAN(der_f *ihs):unary_opr(ihs){}
	der_f* _der(VAR id) const{
		return new _MUL(
			new _INV(
				new _MUL(
					new _COS(ihs->clone()),
					new _COS(ihs->clone())
				)
			),
			ihs->der(id)
		);
	};
	der_f* clone() const{
		return new _TAN(ihs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return tan(ihs->value_at(vmap));
	};
	~_TAN(){};
};
class _COSH;
class _SINH;
class _COSH: public unary_opr {public:
	string op_sym() const {return "cosh";}
	_COSH(der_f *ihs):unary_opr(ihs){}
	der_f* _der(VAR id) const;
	der_f* clone() const{
		return new _COSH(ihs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return cosh(ihs->value_at(vmap));
	};
	~_COSH(){};
};
class _SINH: public unary_opr {public:
	string op_sym() const {return "sinh";}
	_SINH(der_f *ihs):unary_opr(ihs){}
	der_f* _der(VAR id) const;
	der_f* clone() const{
		return new _SINH(ihs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return sinh(ihs->value_at(vmap));
	};
	~_SINH(){};
};
der_f* _COSH::_der(VAR id) const{
	return new _MUL(
		new _SIN(
			ihs->clone()
		),
		ihs->der(id)
	);
};
der_f* _SINH::_der(VAR id) const{
	return new _MUL(
		new _COSH(ihs->clone()),
		ihs->der(id)
	);
};
class _LN: public unary_opr {public:
	string op_sym() const {return "ln";}
	_LN(der_f *ihs):unary_opr(ihs){}
	der_f* _der(VAR id) const{
		return new _MUL(
			new _INV(ihs->clone()),
			ihs->der(id)
		);
	};
	der_f* clone() const{
		return new _LN(ihs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return log(ihs->value_at(vmap));
	};
	~_LN(){};
};
class _TANH: public unary_opr {public:
	string op_sym() const {return "tanh";}
	_TANH(der_f *ihs):unary_opr(ihs){}
	der_f* _der(VAR id) const{
		return new _MUL(
			new _SUB(
				new _CONSTANT(1),
				new _MUL(
					this->clone(),
					this->clone()
				)
			),
			ihs->der(id)
		);
	};
	der_f* clone() const{
		return new _TANH(ihs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return tanh(ihs->value_at(vmap));
	};
	~_TANH(){};
};

class _POW : public binary_opr{ public:
	string op_sym() const{return "^";}
	_POW(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* _der(VAR id) const{
		return new _MUL(
			this->clone(),
			new _ADD(
				new _MUL(
					rhs->der(id),
					new _LN(lhs->clone())
				),
				new _MUL(
					new _DIV(
						rhs->clone(),
						lhs->clone()
					),
					lhs->der(id)
				)
			)
		);
	};
	der_f* clone() const{
		return new _POW(lhs->clone(),rhs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return pow((lhs->value_at(vmap)),(rhs->value_at(vmap)));
	}
};

class _LOG : public binary_opr{ public:
	string op_sym() const{return "^";}
	_LOG(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* _der(VAR id) const{
		auto p=new _DIV(new _LN(rhs->clone()),new _LN(lhs->clone()));
		auto res=p->der(id);
		delete p;
		return res;
	};
	der_f* clone() const{
		return new _LOG(lhs->clone(),rhs->clone());
	};
	double value_at(map<VARID, double> vmap) const{
		return log(lhs->value_at(vmap))/log((rhs->value_at(vmap)));
	}
};

int main(int argc, char* argv[]){
	VAR x("x");
	auto s=MUL(
		ADD(
			CONSTANT(1),
			IDENTITY(x)
		),
		SUB(
			ADD(
				CONSTANT(1),
				IDENTITY(x)
			),
			CONSTANT(1)
		)
	);
	cout<<(string)(*s)<<endl;
	auto s_der=s->der(x);
	cout<<(string)(*s_der)<<endl;
	delete s_der;
	delete s;
	cout<<"hw"<<endl;
}
