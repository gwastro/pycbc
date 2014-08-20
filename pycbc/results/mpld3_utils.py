""" This module provides functionality to extend mpld3 
"""
import mpld3, mpld3.plugins, mpld3.utils

class ClickLink(mpld3.plugins.PluginBase):
    """Plugin for following a link on click"""
    
    JAVASCRIPT = """
    mpld3.register_plugin("clicklink", ClickLink);
    ClickLink.prototype = Object.create(mpld3.Plugin.prototype);
    ClickLink.prototype.constructor = ClickLink;
    ClickLink.prototype.requiredProps = ["id"];
    ClickLink.prototype.defaultProps = {
        links: null
    }
    function ClickLink(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    
    ClickLink.prototype.draw = function(){
        var obj = mpld3.get_element(this.props.id);
        var links = this.props.links;
        
        obj.elements().on("mousedown",
                          function(d, i){ 
                                           window.open(links[i]);
                                        }
                          );
    }
    """
    def __init__(self, points, links):
        self.dict_ = {"type": "clicklink",
                      "id": mpld3.utils.get_id(points),
                      "links": links,
                      }
