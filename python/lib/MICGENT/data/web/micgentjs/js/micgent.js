/*!
 * MICGENT JavaScript Library v1.0
 * http://andreyto.github.com/mgtaxa/
 *
 * Copyright 2018-2020, AstraZeneca
 * Licensed under the GPL Version 3 license.
 *
 */

//namespace as a self-executing function
var micgent = (function () {

    function callName(name, context /*, args */) {
        if(!context) {
            context = window;
        }
        var args = Array.prototype.slice.call(arguments, 2);
        var namespaces = name.split(".");
        var func = namespaces.pop();
        for (var i = 0; i < namespaces.length; i++) {
            context = context[namespaces[i]];
        }
        return context[func].apply(context, args);
    }

    function popoverDataAttributes(x) {
        var attrs = {};
        x.map(function (el) {
            if (el.name) {
                attrs[el.name] = el.value;
            }
        });
        return attrs;
    }

    function initBrowser(cont, config, succeed, fail, context) {
        if(!context) {
            context = window;
        }
        $.getJSON(config, function (config) {
            $("#title_geno_browser").html(config["title"]);
            $("#descr_geno_browser").html(config["descr"]);
            var div = $(cont)[0];
            var igv_opt = config["igv"];
            igv.createBrowser(div, igv_opt)
            .then(function(browser) {
            browser.on('trackclick', function (track, popoverData) {
                track_clicked = track;
                if (track.config.onClick) {
                    return callName(track.config.onClick, context, browser, track, popoverData);
                }
                return undefined;
            });
            if (succeed)
                succeed(browser);
        })}).fail(function (jqxhr, textStatus, error) {
            var err = textStatus + ", " + error;
            console.log("Request Failed: " + err);
            if (fail)
                fail(textStatus, error);
        });
    }

    //' This selects the the row (vizually by default), and
    //' returns data of selected rows [{}].
    //' Return matching row data object or null on non-match,
    //' or undefined value if Tabulator has not been initialized yet on this
    //' element.
    function tabulatorSelectAndGetData(tbl, selector, select = true, scrollTo = true, onActivateQueue = null) {
        var ret = undefined;
        var tbl = $(tbl);
        if (isTabulator(tbl)) {
            //if key did not match anything, getSelectedData will be returning []
            var row = tbl.tabulator("getRow", selector);
            ret = null;
            if (row) {
                ret = row.getData();
                //scroll only works right when visible, otherwise corrupts the view
                if (scrollTo) {
                    var action = function () {
                        row.scrollTo();
                    }
                    if (onActivateQueue) {
                        //onActivateQueue.push(action);
                    } else {
                        action();
                    }
                }
            }
            if (select) {
                tbl.tabulator("selectRow", selector);
            }
        }
        return ret;
    }

    function setDefaults(options, defaultOptions) {
        for (var index in defaultOptions) {
            if (typeof options[index] == "undefined")
                options[index] = defaultOptions[index];
        }
    }

    function isNull(data) {
        return typeof data === 'undefined' || data == null;
    }

    var AsmSampleBrowser = function (cont, configFile) {
        if (cont) {
            cont = $(cont);
        } else {
            cont = $("body");
        }
        this.cont = cont;
        this.configFile = configFile || "browserConf.json";
        var self = this;
        $.getJSON(self.configFile, function (config) {
            config.browserConf = self.configFile; // IGV will open this again
            self.config = config;
            self.initDOM();
            self.initObj();
        }).fail(function (jqxhr, textStatus, error) {
            var err = textStatus + ", " + error;
            console.log("Request Failed: " + err);
        });            
    }

    AsmSampleBrowser.prototype.initDOM = function () {
        var config = this.config;
        var ht = `
                 <div id="geno_browser" style="height: 1000px;min-height: 1000px;"></div>
                 `;
        $(this.cont).append(ht);
    }

    AsmSampleBrowser.prototype.initObj = function () {
        var self = this;
        config = self.config;
        self.initReadsGenoBrowser();
    }

    AsmSampleBrowser.prototype.initReadsGenoBrowser = function () {
        var self = this;
        var contGeno = self.findChild("#geno_browser");
        //contGeno.height($(self.cont).innerHeight());
        if (self.geno_browser === undefined) {

            self.geno_browser = null;

            initBrowser(contGeno, self.config.browserConf,
                    function (browser) {
                        self.geno_browser = browser;
                    },null,self);
        }

    }

    AsmSampleBrowser.prototype.findChild = function (expr) {
        return this.cont.find(expr);
    }

    //-----------------------------------------------
    //--ASMBROWSER START------------------------------------
    //-----------------------------------------------
    var AsmBrowser = function (cont, configFile) {
        if (cont) {
            cont = $(cont);
        } else {
            cont = $("body");
        }
        this.cont = cont;
        this.configFile = configFile || "asmBrowser.json"
        var self = this;

        $.getJSON(self.configFile, function (config) {

            self.config = config

            self.initDOM();
            self.initTabs();
            self.initCommon();
        }).fail(function (jqxhr, textStatus, error) {
            var err = textStatus + ", " + error;
            console.log("Request Failed: " + err);
        });
    }

    AsmBrowser.prototype.initDOM = function () {
        var config = this.config;
        var c = this.cont;
        //multiqc will only show dynamic plots in the iFrame when on CORS server
        //Galaxy is CORS; npm http-server is not.
        //multiqc.html also needs the iFrameResizer JS pilot injected as per docs.
        var ht =
                `<div id="tabs" style="height: 1000px;min-height: 1000px;">
                   <ul>
                       <li id="li_tbl_final"><a href="#tab_tbl_final"><span>Final Assemblies</span></a></li>
                       <li id="li_tbl_all"><a href="#tab_tbl_all"><span>All Assemblies</span></a></li>
                       <li id="li_tbl_sum"><a href="#tab_tbl_sum"><span>Sample Summaries</span></a></li>
                       <li id="li_tbl_dict"><a href="#tab_tbl_dict"><span>Manifest Dictionary</span></a></li>
                       <li id="li_browser"><a href="#tab_geno_browser"><span>Sequences</span></a></li>
                       <li id="li_multiqc"><a href="#tab_multiqc"><span>MultiQC</span></a></li>
                       <li id="li_toc"><a href="toc.html"><span>Downloads</span></a></li>
                   </ul>
                   <div id="tab_tbl_all">
                        <div>Contigs before applying the post-assembly filter.
                        Click on any row for links to genome browser views
                        of the individual reads.
                        Click <a href="#anch_tbl_all_pivot">here</a> to scroll horizontally for more table fields or 
                        to aggregate data in a pivot table.</div>
                        <br/>
                        <div id="tbl_all"></div>
                        <br/>
                        <a name="anch_tbl_all_pivot"></a>
                        <div>Drag fields to re-group the Pivot table</div>                        
                        <div id="tbl_all_pivot"></div>
                   </div>
                   <div id="tab_tbl_final">
                        <div>Contigs that passed the post-assembly filter.
                        Click on any row for links to genome browser views
                        of the individual reads.
                        Click <a href="#anch_tbl_final_pivot">here</a> to scroll horizontally for more table fields or 
                        to aggregate data in a pivot table.</div>
                        <br/>
                        <div id="tbl_final"></div>
                        <br/>
                        <a name="anch_tbl_final_pivot"></a>
                        <div>Drag fields to re-group the Pivot table</div>
                        <div id="tbl_final_pivot"></div>
                    </div>
                   <div id="tab_tbl_sum">
                        <div>Per-sample summary table.
                        Click <a href="#anch_tbl_sum_pivot">here</a> to scroll horizontally for more table fields or 
                        to aggregate data in a pivot table.</div>
                        <br/>
                        <div id="tbl_sum"></div>
                        <br/>
                        <a name="anch_tbl_sum_pivot"></a>
                        <div>Drag fields to re-group the Pivot table</div>
                        <div id="tbl_sum_pivot"></div>
                    </div>                    
                   <div id="tab_tbl_dict">
                        <div>Description of the fields from the Assembly and Summary tables</div>
                        <br/>
                        <div id="tbl_dict"></div>
                        <br/>
                    </div>                                        
                   <div id="tab_geno_browser">
                        <div>Contigs that passed the post-assembly filter are shown as aligned against the
                        provided common reference(s). Each contig is shown as one "read". Note: this genome
                        browser does not highlight ambiguous bases in the contigs as it does for standard
                        ACTG nucleotides.
                        Click on any contig for links to genome browser views
                        of the individual sequencing reads.                        
                        </div>
                        <br/>
                        <div id="geno_browser"></div>
                   </div>
                   <div id="tab_multiqc">
                        <style>iframe{width: 1px;min-width: 100%;}</style>
                        <iframe id="multiqc" src="${config.multiqc.html}" scrolling="no"></iframe>
                        <script>iFrameResize({checkOrigin:false}, '#multiqc')</script>        
                   </div>
                 </div>
                 <div id="msg_dlg"></div>
                `;
        c.append(ht);
    }

    function replayOnActivateQueue(queue) {
        var any = false;
        var action;
        if (queue) {
            while (action = queue.shift()) {
                any = true;
                action();
            }
        }
        return any;
    }

    AsmBrowser.prototype.initTabs = function () {

        var self = this;

        var config = self.config;
        // Tabulator tables can be only correctly initialized when their 
        // DOM element is fully visible. This causes the convoluted
        // construction process below.
        // self.onActivateQueues = {tbl_all: [], tbl_final: [], geno_browser: []};
        $('#tabs').tabs({
            heightStyle: "fill",
            event: 'click',
            create: function (event, ui) {
                var id_el = $(ui.panel).prop("id").replace("tab_", "");
                if (id_el === "tbl_all") {
                    self.initTbl('tbl_all', config.manifestAll);
                } else if (id_el === "tbl_final") {
                    self.initTbl('tbl_final', config.manifestFinal);
                } else if (id_el === "tbl_sum") {
                    self.initTbl('tbl_sum', config.manifestSum);
                } else if (id_el === "tbl_dict") {
                    self.initTbl('tbl_dict', config.manifestDict);                                        
                } else if (id_el === "geno_browser") {
                    self.initContigGenoBrowser();
                }
            },
            activate: function (event, ui) {
                var id_el = $(ui.newPanel).prop("id").replace("tab_", "");
                var el = self.findChild("#" + id_el);
                if (isTabulator(el)) {
                    //replayOnActivateQueue is useless for now and the queues
                    //should not be populated becaue there is no callback
                    //for the end of redraw(false). Calling replay here corrupts
                    //the table view
                    el.tabulator('redraw', false);
                } else {
                    if (id_el === "tbl_all") {
                        self.initTbl('tbl_all', config.manifestAll);
                    } else if (id_el === "tbl_final") {
                        self.initTbl('tbl_final', config.manifestFinal);
                    } else if (id_el === "tbl_sum") {
                        self.initTbl('tbl_sum', config.manifestSum);
                    } else if (id_el === "tbl_dict") {
                        self.initTbl('tbl_dict', config.manifestDict);                        
                    } else if (id_el === "geno_browser") {
                        self.initContigGenoBrowser();
                    } else if (id_el === "multiqc") {
                        // refresh iframe on activate, otherwise dynamic
                        // plots do not show up
                        $('#'+id_el)[0].src += '';
                    }
                }

            },
            beforeActivate: function (event, ui) {
                self.findChild("#msg_dlg").dialog("close");
            }            
        });
    }

    AsmBrowser.prototype.initCommon = function () {
        var self = this;
        self.findChild("#msg_dlg").dialog({ autoOpen: false });
    }

    function isTabulator(selector) {
        return $(selector).attr("class") === "tabulator";
    }

    AsmBrowser.prototype.initContigGenoBrowser = function () {
        var self = this;
        var contGeno = self.findChild("#geno_browser");
        //contGeno.height($(self.cont).innerHeight());
        if (self.geno_browser === undefined) {

            self.geno_browser = null;

            initBrowser(contGeno, self.config.browserConf,
                    function (browser) {
                        self.geno_browser = browser;
                    },null,self);
        }

    }

    AsmBrowser.prototype.initTbl = function (name, dataFile) {
        var self = this;

        $.get(dataFile, function (input) {
           //alternative CSV parser: https://www.papaparse.com/
           var data = $.csv.toObjects(input, {"separator": "\t"});
            var columns = [];
            if (data.length>0) {
                var fields = Object.keys(data[0]);
                columns = fields.map(function (x) {
                    var rec = {title: x, field: x};
                    if (x === "SampleID" || x === "SeqID" || x === "Asm_Failed") {
                        rec.frozen = true;
                        rec.headerFilter = "input";
                    }
                    if(["file1","file2","Asm_sig_inp1","Asm_sig_inp2"].includes(x)) {
                        rec.width = 100;
                    }
                    return rec;
                });
            }
            //var onActivateQueue = self.onActivateQueues[name];
            var dom_el = self.findChild("#" + name);
            var index = "SeqID";
            var rowClick = function (e, row) { //trigger an alert message when the row is clicked
                    var SeqID = row.getData().SeqID;
                    var seqData = self.selectSeqIDInTables(SeqID);
                    var msg = self.makeSamplePopupHTML(SeqID, seqData, name);
                    var dlg_el = self.findChild("#msg_dlg");
                    dlg_el.dialog( "option", "title", "Assembly Details" );
                    //var msg = `Row SampleID ${row.getData().SampleID} SeqID ${row.getData().SeqID}`
                    dlg_el.html(msg);
                    dlg_el.dialog("open");
                } 
            var add_pivot = true;
            var pivot_cols = ["Asm_ref_name"];
            var pivot_vals = ["Asm_reads"];
            if(name == "tbl_sum") {
                index = "SampleID"
                rowClick = null;
                pivot_cols = ["Asm_Failed"];
                pivot_vals = ["Asm_SeqCountFinal"];
            } else if(name == "tbl_dict") {
                index = "Variable"
                add_pivot = false;
                rowClick = null;
            }
            dom_el.tabulator({
                height: 800, //A value is needed for speed-up with virtual DOM
                //layout: "fitDataFill",
                //responsiveLayout: true,
                placeholder: "Loading data",
                index: index,
                selectable: 1,
                columns: columns,
                rowClick: rowClick
                //renderComplete: function () {
                //    replayOnActivateQueue(onActivateQueue);
                //}
            });
            dom_el.tabulator("setData", data);
            if(add_pivot) {
                var dom_pivot = self.findChild("#" + name + "_pivot")
                dom_pivot.pivotUI(
                  data, {
                    rows: ["SampleID"],
                    cols: pivot_cols,
                    vals: pivot_vals,
                    aggregatorName: "Average",
                    rendererName: "Heatmap",
                    renderers: $.extend(
                      $.pivotUtilities.renderers,
                      $.pivotUtilities.d3_renderers,
                      $.pivotUtilities.c3_renderers,
                      $.pivotUtilities.plotly_renderers,
                      $.pivotUtilities.export_renderers
                    )
                  });
            }

        }).fail(function (jqxhr, textStatus, error) {
            var err = textStatus + ", " + error;
            console.log("Request Failed: " + err);
        });

    }

    AsmBrowser.prototype.findChild = function (expr) {
        return this.cont.find(expr);
    }

    AsmBrowser.prototype.makeSamplePopupHTML = function (SeqID, seqData, tbl_name="tbl_final") {
        var self = this;
        var config = self.config;
        //should by default pull data from whatever table is the first tab - otherwise
        //it will not be available until user activates the tab
        var SampleID = seqData[tbl_name].SampleID;
        var Asm_SeqStatus = seqData[tbl_name].Asm_SeqStatus;
        var base_msg = `
            Sample ID: ${SampleID}<br/>        
            `;
        var asm_msg = "Assembly process failed to produce sequences";
        if(Asm_SeqStatus !== "miss") {
            var cluster = seqData[tbl_name].Asm_cluster;
            asm_msg = `
                Contig ID: ${SeqID}<br/>
                <a href="${config.samplesDir}/${SampleID}/clusters/${cluster}/asm_map.html" rel="noopener noreferrer" target="_blank">
                New Window with Reads Mapped to Assembly</a><br/>
                <a href="${config.samplesDir}/${SampleID}/clusters/${cluster}/ref_map.html" rel="noopener noreferrer" target="_blank">
                New Window with Reads Mapped to Reference</a>
                `;
        }
        return '<div>' + base_msg + asm_msg + '</div>';
    }

    AsmBrowser.prototype.selectSeqIDInTables = function (SeqID) {
        var self = this;
        var ret = {};
        for (var tbl in {"tbl_all": null, "tbl_final": null}) {
            ret[tbl] = tabulatorSelectAndGetData(self.findChild("#" + tbl), SeqID);
        }
        return ret;
    }

    AsmBrowser.prototype.contigsOnClick = function (browser, track, popoverData) {
        var self = this;
        var attrs = popoverDataAttributes(popoverData);
        var ret = undefined;
        var idField = "Read Name";
        var SeqID = attrs[idField];
        if (SeqID) { // did not click coverage hist or empty space
            var config = self.config;
            var seqData = self.selectSeqIDInTables(SeqID);
            seqData["seq"] = attrs;
            ret = self.makeSamplePopupHTML(SeqID, seqData);
        }

        return ret; // undefined keeps default pop-over
    }


    return {
        AsmBrowser: AsmBrowser,
        AsmSampleBrowser: AsmSampleBrowser
    }
    //	end of namespace definition
})();
