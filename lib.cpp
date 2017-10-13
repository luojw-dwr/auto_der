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
	static map<int, VAR> vars;
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
		return (vid==rhs.vid)&&(vname==rhs.vname);
	};
	bool operator<(VAR rhs) const{
		return (vname==rhs.vname)&&(vid==rhs.vid);
	}
};
int VAR::vid_acc=0;

class der_f{ public:
	virtual operator string() const =0;
	virtual long double value_at(map<VARID, long double>) const =0;
	long double value_at_origin() const{
		map<VARID, long double> empty_vmap;
		return value_at(empty_vmap);
	}
	virtual der_f* der(VAR) const =0;
	virtual der_f* clone() const =0;
	virtual ~der_f(){};
	virtual bool has_var(VAR) const =0;
	virtual bool is_constant() const {return false;}//Not so precise
};

class CONSTANT : public der_f { public:
	long double v;
	CONSTANT(long double v):v(v){};
	operator string() const {
		return std::to_string(v);
	}
	der_f* der(VAR id) const{
		return new CONSTANT(0);
	};
	der_f* clone() const{
		return new CONSTANT(v);
	};
	bool has_var(VAR id) const{return false;};
	long double value_at(map<VARID, long double> vmap) const{
		return v;
	}
	bool is_constant() const{
		return true;
	}
};

class IDENTITY : public der_f { public:
	VAR id;
	IDENTITY(VAR id){
		id.clone(&this->id);
	};
	operator string() const {
		return id.vname;
	}
	der_f* der(VAR _id) const{
		if(id.eq(_id)){
			return new CONSTANT(1);
		} else {
			return new CONSTANT(0);
		}
	};
	der_f* clone() const{
		return new IDENTITY(id);
	};
	
	bool has_var(VAR _id) const{return id.eq(_id);};
	long double value_at(map<VARID, long double> vmap) const{
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
			this->ihs=new CONSTANT(ihs->value_at_origin());
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
			this->lhs=new CONSTANT(lhs->value_at_origin());
			delete lhs;
		}else{
			this->lhs=lhs;
		}
		if(rhs->is_constant()){
			this->rhs=new CONSTANT(rhs->value_at_origin());
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
class ADD : public binary_opr{ public:
	string op_sym() const{return "+";}
	ADD(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* der(VAR id) const{
		return new ADD(lhs->der(id),rhs->der(id));
	};
	der_f* clone() const{
		return new ADD(lhs->clone(),rhs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return (lhs->value_at(vmap))+(rhs->value_at(vmap));
	}
};
class SUB : public binary_opr{ public:
	string op_sym() const{return "-";}
	SUB(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* der(VAR id) const{
		return new SUB(lhs->der(id),rhs->der(id));
	};
	der_f* clone() const{
		return new SUB(lhs->clone(),rhs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return (lhs->value_at(vmap))-(rhs->value_at(vmap));
	}
};
class MUL : public binary_opr{ public:
	string op_sym() const{return "*";}
	MUL(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* der(VAR id) const{
		return new ADD(
			new MUL(lhs->der(id),rhs->clone()),
			new MUL(lhs->clone(),rhs->der(id))
		);
	};
	der_f* clone() const{
		return new MUL(lhs->clone(),rhs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return (lhs->value_at(vmap))*(rhs->value_at(vmap));
	}
};
class DIV : public binary_opr{ public:
	string op_sym() const{return "/";}
	DIV(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* der(VAR id) const{
		return new DIV(
			new SUB(
				new MUL(lhs->der(id),rhs->clone()),
				new MUL(lhs->clone(),rhs->der(id))
			),
			new MUL(
				rhs->clone(),
				rhs->clone()
			)
		);
	};
	der_f* clone() const{
		return new DIV(lhs->clone(),rhs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return (lhs->value_at(vmap))/(rhs->value_at(vmap));
	}
};

class NEG: public unary_opr {public:
	string op_sym() const {return "-";}
	NEG(der_f *ihs):unary_opr(ihs){}
	der_f* der(VAR id) const{
		return new NEG(ihs->der(id));
	}
	der_f* clone() const{
		return new NEG(ihs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return -(ihs->value_at(vmap));
	}
};
class INV: public unary_opr {public:
	string op_sym() const {return "inv";}
	INV(der_f *ihs):unary_opr(ihs){}
	der_f* der(VAR id) const{
		return new NEG(
			new DIV(
				ihs->der(id),
				new MUL(
					ihs->clone(),
					ihs->clone()
				)
			)
		);
	}
	der_f* clone() const{
		return new INV(ihs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return 1.0/(ihs->value_at(vmap));
	}
};
class EXP: public unary_opr {public:
	string op_sym() const {return "exp";}
	EXP(der_f *ihs):unary_opr(ihs){}
	der_f* der(VAR id) const{
		return new MUL(
			this->clone(),
			ihs->der(id)
		);
	}
	der_f* clone() const{
		return new EXP(ihs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return exp(ihs->value_at(vmap));
	};
	~EXP(){};
};

class COS: public unary_opr {public:
	string op_sym() const {return "cos";}
	COS(der_f *ihs):unary_opr(ihs){}
	der_f* der(VAR id) const{
		return new MUL(
			new NEG(
				new COS(
					new SUB(
						new CONSTANT(M_PI/2),
						ihs->clone()
					)
				)
			),
			ihs->der(id)
		);
	}
	der_f* clone() const{
		return new COS(ihs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return cos(ihs->value_at(vmap));
	};
	~COS(){};
};
class SIN: public unary_opr {public:
	string op_sym() const {return "sin";}
	SIN(der_f *ihs):unary_opr(ihs){}
	der_f* der(VAR id) const{
		return new MUL(
			new COS(ihs->clone()),
			ihs->der(id)
		);
	};
	der_f* clone() const{
		return new SIN(ihs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return sin(ihs->value_at(vmap));
	};
	~SIN(){};
};
class LN: public unary_opr {public:
	string op_sym() const {return "ln";}
	LN(der_f *ihs):unary_opr(ihs){}
	der_f* der(VAR id) const{
		return new MUL(
			new INV(ihs->clone()),
			ihs->der(id)
		);
	};
	der_f* clone() const{
		return new LN(ihs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return log(ihs->value_at(vmap));
	};
	~LN(){};
};

class POW : public binary_opr{ public:
	string op_sym() const{return "^";}
	POW(der_f* lhs, der_f* rhs):binary_opr(lhs,rhs){};
	der_f* der(VAR id) const{
		return new MUL(
			this->clone(),
			new ADD(
				new MUL(
					rhs->der(id),
					new LN(lhs->clone())
				),
				new MUL(
					new DIV(
						rhs->clone(),
						lhs->clone()
					),
					lhs->der(id)
				)
			)
		);
	};
	der_f* clone() const{
		return new POW(lhs->clone(),rhs->clone());
	};
	long double value_at(map<VARID, long double> vmap) const{
		return pow((lhs->value_at(vmap)),(rhs->value_at(vmap)));
	}
};

int main(int argc, char* argv[]){
	VAR x("x");
	auto s=new POW(new LN(new ADD(new IDENTITY(x),new CONSTANT(5))),new EXP(new SIN(new NEG(new IDENTITY(x)))));
	cout<<s->value_at_origin()<<endl;
	map<VARID, long double> vmap;
	vmap[x.vid]=3;
	cout<<s->value_at(vmap)<<endl;
	cout<<(string)(*s)<<endl;
	cout<<s->value_at_origin()<<endl;
	auto s_der=s->der(x);
	cout<<(string)(*s_der)<<endl;
	delete s_der;
	delete s;
	cout<<"hw"<<endl;
}
