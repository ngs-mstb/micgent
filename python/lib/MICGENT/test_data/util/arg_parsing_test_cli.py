from MICGENT import yaml_util
import argh

@argh.arg("--config",type=yaml_util.load_yaml)
def my_g(x,y="a",config={},out="my_g_out.yaml"):
    ret = dict(x=x,y=y,config=config,out=out)
    yaml_util.dump_yaml(ret,out)
    return ret

## import package module and add argh entry points

def _main():
    from MICGENT import arg_parsing
    parser = arg_parsing.ArghParserChainedConfig()
    parser.add_commands([
        my_g
    ])
    parser.dispatch()


if __name__ == "__main__":
    _main()
